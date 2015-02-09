ocamlopt -annot -o main.opt -I `opam config var lib`/sundialsml bigarray.cmxa sundials.cmxa graphics.cmxa base_model.ml
# ocamlopt -inline 100 -annot -o main.opt -I `opam config var lib`/sundialsml bigarray.cmxa sundials.cmxa graphics.cmxa diffusion.ml
