import Module from './coolprop.js'
var coolprop = await Module();
console.log('32 F in K:', coolprop.F2K(32.0)); // values can be returned, etc.
console.log('NBP of water in K:', coolprop.PropsSI('T','P',101325,'Q',0,'Water')); // values can be returned, etc.

var Tnbp = coolprop.PropsSI('T','P',101325,'Q',0,'Water');
if (Tnbp < 373.124 || Tnbp > 373.125){
  console.error("NBP is invalid");
  process.exit(1);
}