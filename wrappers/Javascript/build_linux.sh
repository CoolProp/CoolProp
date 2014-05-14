em++ -O2 ../../CoolProp/*.cpp -I../../CoolProp -c  -s EXPORTED_FUNCTIONS="['_main','_F2K','_Props1','_PropsS']"
em++ -O2 -o coolprop.js ../../CoolProp/CoolPropDLL.cpp *.o -I../../CoolProp -DEXTERNC  -s EXPORTED_FUNCTIONS="['_main','_F2K','_Props1','_PropsS']"

# Can only use compression with HTML :(
# --compression ~/Code/emscripten/third_party/lzma.js/lzma-native,/home/xubuntu/Code/emscripten/third_party/lzma.js/lzma-decoder.js,LZMA.decompress

# Using closure compiler to compress javascript file
#java -jar compiler.jar --js precursor.js --js_output_file coolprop.js --compilation_level ADVANCED_OPTIMIZATIONS --language_in ECMASCRIPT5

rm *.o