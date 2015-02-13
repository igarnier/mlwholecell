(* We need a minimum of structure on species. *)

module type Species =
  sig

    type t

    val compare   : t -> t -> int
    val to_string : t -> string

  end

module MakeOde(S : Species) =
  struct

    type multiset = (S.t * int) list

    type reaction =
      { input  : multiset;
        output : multiset;
        rate   : ((S.t -> float) -> float)
      }

    type pre_reaction = multiset * multiset

    (* let rec tree_to_multiset = *)
    (*   function *)
    (*   | Join(lhs, rhs) -> *)
    (*     (tree_to_multiset lhs) @ (tree_to_multiset rhs) *)
    (*   | Data(s,n) ->  *)
    (*     [(s,n)] *)

    (* Make sure the multiset is well-formed, i.e. no redundancy. *)
    let rec normalize_multiset_aux l =
      match l with
      | [] | [_] -> l
      | (s1, n1) :: (s2, n2) :: tl ->
        if S.compare s1 s2 = 0 then
          normalize_multiset_aux ((s1, n1 + n2) :: tl)
        else
          (s1,n1) :: (normalize_multiset_aux ((s2, n2) :: tl))

    let normalize_multiset l =
      let sorted = List.sort (fun (s1, _) (s2, _) -> S.compare s1 s2) l in
      normalize_multiset_aux sorted

    let print_multiset l = 
      List.fold_left (fun acc (s, n) ->
        Printf.sprintf "%s (%s, %d)" acc (S.to_string s) n
      ) "" l

    (* TODO: not so easy to print rates without some kind of metaprogramming. *)
    let print_reaction r =
      Printf.sprintf "%s --> %s" (print_multiset r.input) (print_multiset r.output)

    (* Lighteight syntax for multisets reactions *)

    let (&) : multiset -> S.t * int -> multiset = fun tl sn -> sn :: tl

    let (==>) lhs rhs = (lhs, rhs)
    let (<==) lhs rhs = (rhs, lhs)

    let (@@) (lhs, rhs) rate =
      { input  = lhs;
        output = rhs;
        rate   = rate
      }

    let mset x = [x]

  (* From a bunch of reactions, one easily produces a matrix shaped field *)
    let to_ode reactions = TODO

  end

module Species =
struct 
  
  type t = A | B | C 

  let compare : t -> t -> int = compare

  let to_string = function
    | A -> "A"
    | B -> "B"
    | C -> "C"

end

module R = MakeOde(Species)

open Species
open R

let reac  = ((mset (A, 1) & (B, 3)) ==> (mset (A, 2) & (B, 2))) @@ (fun k -> (k A) +. (k B))
let ireac = ((mset (A, 1) & (B, 3)) <== (mset (A, 2) & (B, 2))) @@ (fun k -> (k A) +. (k B))
 
