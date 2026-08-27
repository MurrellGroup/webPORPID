declare module "@biowasm/aioli" {
  export default class Aioli {
    constructor(tools: Array<string | { tool: string; program: string; version: string; urlPrefix: string }>);
    write(options: { path: string; buffer: Uint8Array }): Promise<void>;
    exec(command: string): Promise<string>;
  }
}
