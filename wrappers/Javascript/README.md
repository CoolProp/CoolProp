# CoolProp

Official WebAssembly build of [CoolProp](https://coolprop.github.io/), an
open-source database of thermophysical properties for pure fluids,
pseudo-pure fluids, mixtures and humid air.

The package ships the emscripten build as an ES6 module (`coolprop.js`) plus
its WebAssembly binary (`coolprop.wasm`), usable in Node.js, bundlers
(Vite, webpack, ...) and the browser.

## Install

```bash
npm install coolprop
```

## Usage

```js
import Module from 'coolprop';

// The default export is the emscripten module factory
const coolprop = await Module();

// High-level interface
console.log(coolprop.F2K(32));                                   // 273.15
console.log(coolprop.PropsSI('T', 'P', 101325, 'Q', 0, 'Water')); // 373.124...

// Low-level AbstractState interface
const AS = coolprop.factory('HEOS', 'Water');
AS.update(coolprop.input_pairs.PQ_INPUTS, 101325, 0);
console.log('T:', AS.T(), 'rho:', AS.rhomass());
AS.delete();
```

## Bundlers and asset hosting

The module locates `coolprop.wasm` relative to itself with
`import.meta.url`, which Vite and webpack understand natively. If you serve
the assets from a custom location, override the resolution with
`locateFile`:

```js
const coolprop = await Module({ locateFile: (file) => `/wasm/${file}` });
```

When self-hosting, make sure the `.wasm` extension is served with the MIME
type `application/wasm`.

## API

Full documentation of the high-level (`PropsSI`, `HAPropsSI`, ...) and
low-level (`AbstractState`) interfaces is at
<https://coolprop.github.io/coolprop/>, and a live demo in the
[CoolPropJavascriptDemo](https://github.com/dvd101x/CoolPropJavascriptDemo)
repository. TypeScript declarations are included (`coolprop.d.ts`).

## License

MIT -- see the LICENSE file. Source code and issue tracker:
<https://github.com/CoolProp/CoolProp>
