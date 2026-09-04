type WorkerRequest = Record<string, unknown> & { type: string };

class WorkerClient {
  private nextId = 1;
  private pending = new Map<number, { resolve(value: unknown): void; reject(error: Error): void }>();
  private tail: Promise<unknown> = Promise.resolve();
  private failure?: Error;
  constructor(readonly worker: Worker) {
    worker.addEventListener("message", (event: MessageEvent<{ id: number; result?: unknown; error?: string }>) => {
      const pending = this.pending.get(event.data.id); if (!pending) return; this.pending.delete(event.data.id);
      if (event.data.error) pending.reject(new Error(event.data.error)); else pending.resolve(event.data.result);
    });
    worker.addEventListener("error", (event) => {
      this.failure = new Error(event.message || "A compute worker crashed.");
      for (const pending of this.pending.values()) pending.reject(this.failure);
      this.pending.clear();
    });
  }

  call<T>(request: WorkerRequest, transfer: Transferable[] = []): Promise<T> {
    if (this.failure) return Promise.reject(this.failure);
    const task = () => new Promise<T>((resolve, reject) => {
      const id = this.nextId++; this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject });
      this.worker.postMessage({ id, ...request }, { transfer });
    });
    const result = this.tail.then(task, task); this.tail = result.catch(() => undefined); return result;
  }
  async terminate() {
    try { await this.call({ type: "shutdown" }); }
    catch { /* a crashed worker has no heap left to release cooperatively */ }
    finally { this.worker.terminate(); }
  }
}

export class CoreWorkerPool {
  private cursor = 0;
  private closed = false;
  private constructor(readonly clients: WorkerClient[]) {}

  static async create(size: number, module: WebAssembly.Module, config: Uint8Array) {
    const clients = Array.from({ length: Math.max(1, size) }, () => new WorkerClient(new Worker(new URL("./core-worker.ts", import.meta.url), { type: "module" })));
    try {
      await Promise.all(clients.map((client) => {
        const copy = config.slice().buffer; return client.call({ type: "initialize", module, config: copy }, [copy]);
      }));
    } catch (cause) { await Promise.all(clients.map((client) => client.terminate())); throw cause; }
    return new CoreWorkerPool(clients);
  }

  any<T>(request: WorkerRequest, transfer: Transferable[] = []) {
    const client = this.clients[this.cursor++ % this.clients.length]; return client.call<T>(request, transfer);
  }
  at<T>(index: number, request: WorkerRequest, transfer: Transferable[] = []) { return this.clients[index].call<T>(request, transfer); }
  async close() { if (this.closed) return; this.closed = true; await Promise.all(this.clients.map((client) => client.terminate())); }
}
