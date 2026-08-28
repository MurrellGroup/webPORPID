# Third-party notices

This project vendors or packages the following third-party components:

- Alivibe-compatible MSA WebAssembly source, the pinned Alivibe editor/round-trip bridge, and Newick layout/rendering code from the MIT-licensed MurrellGroup/swig repository, vendored as standalone code under `wasm/vendor/`, `public/tools/`, and `src/`.
- FastTree WebAssembly assets from the BioWASM/Aioli ecosystem, packaged under `public/biowasm/fasttree/` and `cli/assets/`.
- React, Vite, TypeScript, fflate, MessagePack, and other npm dependencies listed in `package.json`.

The generated release zip excludes `node_modules`, build caches, upstream validation clones, and quarantined original-demo material.

If publishing compiled FastTree artifacts, include the matching upstream FastTree and BioWASM license notices required by the artifact source used for release.
