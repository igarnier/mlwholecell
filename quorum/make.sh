ocamlopt -w A -annot -o main.opt -I `opam config var lib`/sundialsml bigarray.cmxa sundials.cmxa graphics.cmxa main.ml
# ocamlopt -inline 100 -annot -o main.opt -I `opam config var lib`/sundialsml bigarray.cmxa sundials.cmxa graphics.cmxa diffusion.ml
