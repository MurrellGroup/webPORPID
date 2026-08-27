type WorkerRequest = Record<string, unknown> & { type: string };

class WorkerClient {
  private nextId = 1;
  private pending = new Map<number, { resolve(value: unknown): void; reject(error: Error): void }>();
  private tail: Promise<unknown> = Promise.resolve();
  constructor(readonly worker: Worker) {
    worker.addEventListener("message", (event: MessageEvent<{ id: number; result?: unknown; error?: string }>) => {
      const pending = this.pending.get(event.data.id); if (!pending) return; this.pending.delete(event.data.id);
      if (event.data.error) pending.reject(new Error(event.data.error)); else pending.resolve(event.data.result);
    });
    worker.addEventListener("error", (event) => {
      for (const pending of this.pending.values()) pending.reject(new Error(event.message || "A compute worker crashed."));
      this.pending.clear();
    });
  }

  call<T>(request: WorkerRequest, transfer: Transferable[] = []): Promise<T> {
    const task = () => new Promise<T>((resolve, reject) => {
      const id = this.nextId++; this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject });
      this.worker.postMessage({ id, ...request }, { transfer });
    });
    const result = this.tail.then(task, task); this.tail = result.catch(() => undefined); return result;
  }
  terminate() { this.worker.terminate(); }
}

export class CoreWorkerPool {
  private cursor = 0;
  private constructor(readonly clients: WorkerClient[]) {}

  static async create(size: number, module: WebAssembly.Module, config: Uint8Array) {
    const clients = Array.from({ length: Math.max(1, size) }, () => new WorkerClient(new Worker(new URL("./core-worker.ts", import.meta.url), { type: "module" })));
    await Promise.all(clients.map((client) => {
      const copy = config.slice().buffer; return client.call({ type: "initialize", module, config: copy }, [copy]);
    }));
    return new CoreWorkerPool(clients);
  }

  any<T>(request: WorkerRequest, transfer: Transferable[] = []) {
    const client = this.clients[this.cursor++ % this.clients.length]; return client.call<T>(request, transfer);
  }
  at<T>(index: number, request: WorkerRequest, transfer: Transferable[] = []) { return this.clients[index].call<T>(request, transfer); }
  close() { this.clients.forEach((client) => client.terminate()); }
}
