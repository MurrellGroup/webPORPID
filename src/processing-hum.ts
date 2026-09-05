/**
 * A deliberately audible, very quiet processing hum.
 *
 * Chrome gives recently audible pages more permissive background timer
 * scheduling. A zero-valued or silent track does not qualify, so this graph
 * emits a real signal. It is an optional best-effort aid only: browsers and
 * operating systems remain free to freeze or discard a page.
 */

export const PROCESSING_HUM_LEVEL = 0.007;
const FADE_SECONDS = 0.12;

type AudioContextConstructor = new (options?: AudioContextOptions) => AudioContext;
type AudioGlobal = typeof globalThis & { webkitAudioContext?: AudioContextConstructor };

function browserAudioContext(): AudioContext {
  const audioGlobal = globalThis as AudioGlobal;
  const Constructor = globalThis.AudioContext ?? audioGlobal.webkitAudioContext;
  if (!Constructor) throw new Error("Web Audio is not available in this browser.");
  return new Constructor({ latencyHint: "playback" });
}

interface HumSession {
  context: AudioContext;
  output: GainNode;
  sources: OscillatorNode[];
}

export class ProcessingHum {
  private session?: HumSession;
  private readonly createContext: () => AudioContext;

  constructor(createContext: () => AudioContext = browserAudioContext) { this.createContext = createContext; }

  get active() { return Boolean(this.session); }

  async start(): Promise<void> {
    if (this.session) {
      if (this.session.context.state === "suspended") await this.session.context.resume();
      return;
    }

    const context = this.createContext();
    try {
      const now = context.currentTime;
      const output = context.createGain();
      output.gain.setValueAtTime(0, now);
      output.gain.linearRampToValueAtTime(PROCESSING_HUM_LEVEL, now + FADE_SECONDS);

      const lowPass = context.createBiquadFilter();
      lowPass.type = "lowpass";
      lowPass.frequency.setValueAtTime(285, now);
      lowPass.Q.setValueAtTime(0.7, now);
      lowPass.connect(output);
      output.connect(context.destination);

      const makeTone = (type: OscillatorType, frequency: number, level: number) => {
        const source = context.createOscillator(), gain = context.createGain();
        source.type = type;
        source.frequency.setValueAtTime(frequency, now);
        gain.gain.setValueAtTime(level, now);
        source.connect(gain); gain.connect(lowPass); source.start(now);
        return source;
      };

      // Low fundamentals plus a tiny slow modulation make a gentle engine-like
      // sound while keeping the master signal around -43 dBFS.
      const sources = [makeTone("triangle", 82, 0.65), makeTone("sine", 123, 0.2)];
      const lfo = context.createOscillator(), lfoDepth = context.createGain();
      lfo.type = "sine"; lfo.frequency.setValueAtTime(0.32, now);
      lfoDepth.gain.setValueAtTime(0.0012, now);
      lfo.connect(lfoDepth); lfoDepth.connect(output.gain); lfo.start(now); sources.push(lfo);

      this.session = { context, output, sources };
      await context.resume();
      if (this.session?.context === context && context.state !== "running")
        throw new Error("The browser did not allow the processing sound to start.");
    } catch (cause) {
      if (this.session?.context === context) this.session = undefined;
      try { await context.close(); } catch { /* best-effort cleanup */ }
      throw cause;
    }
  }

  async stop(): Promise<void> {
    const session = this.session;
    this.session = undefined;
    if (!session) return;
    const { context, output, sources } = session, now = context.currentTime;
    try {
      output.gain.cancelScheduledValues(now);
      output.gain.setValueAtTime(Math.max(0, output.gain.value), now);
      output.gain.linearRampToValueAtTime(0, now + FADE_SECONDS);
      for (const source of sources) try { source.stop(now + FADE_SECONDS); } catch { /* already stopped */ }
      await new Promise<void>((resolve) => globalThis.setTimeout(resolve, Math.ceil((FADE_SECONDS + 0.03) * 1000)));
    } finally {
      try { await context.close(); } catch { /* best-effort cleanup */ }
    }
  }
}
