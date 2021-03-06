module RealArray = Sundials.RealArray
(* module Roots = Sundials.Roots *)

type params = {
  s         : float; (* nutrient concentration *)
  dm        : float; (* mRNA degradation rate *)
  ns        : float; (* nutrient efficiency *)
  nr        : float;
  nt        : float;
  nm        : float;
  nq        : float;
  gamma_max : float; (* max transl. elongation rate *)
  kK_gamma  : float; (* transl. elongation threshold *)
  vt        : float; (* max nutrient import rate *)
  kKt       : float; (* nutient import threshold *)
  vm        : float; (* max enzymatic rate *)
  kKm       : float; (* enzymatic threshold *)
  wr        : float; (* max prot. transcription rate *)
  wt        : float;
  wm        : float;
  wq        : float;
  theta_r   : float; (* ribosome transcription threshold  *)
  theta_nr  : float;
  kKq       : float; (* q-qutoinhibition threshold *)
  hq        : float; (* q-qutoinhibition Hill coef. *)
  k_b       : float; (* mRNA-ribosome binding rate *)
  k_u       : float; (* mRNA-ribosome unbinding rate *)
  mM        : float; (* target mass *)
}

type protein = {
  (* Concentration of prot., free mrna and bound mrna *)
  prot        : float;
  free_mrna   : float;
  bound_mrna  : float;

  (* Rates associated to prot. *)
  nu          : float;
  omega       : float;
}
type cell = {
  s_i    : float;
  a      : float;
  r      : protein;
  e_t    : protein;
  e_m    : protein;
  q      : protein;
}

(* Unused (for now) *)
type ddt_protein = {
  ddt_prot        : float;
  ddt_free_mrna   : float;
  ddt_bound_mrna  : float;
}

(********************************************************)
(* Printing functions *)
(********************************************************)

open Printf

let print x =
  for i = 0 to 13 do
    printf "%f%! " x.{i}
  done;
  printf "\n"

let fprint file =
  let fd = open_out file in
  fun obs x ->
    match x with
    | None -> close_out fd
    | Some data -> fprintf fd "%f\n" (obs data)

let fprint_vec file =
  let fd = open_out file in
  let f (data : Sundials.RealArray.t) =
      (* let bigarr = Nvector.unwrap data in *)
    for i = 0 to Bigarray.Array1.dim data - 1 do
      fprintf fd "%f " data.{i}
    done;
    fprintf fd "\n"
  in
  (fd, f)

(********************************************************)
(* Helper functions *)
(********************************************************)

(* There is room for code factorisation here. *)

let hill vmax kKm x h =
  vmax *. (x ** h) /. (kKm +. (x ** h))

let omega_r a k =
  k.wr *. a /. (a +. k.theta_r)

let omega_m a k =
  k.wm *. a /. (a +. k.theta_nr)

let omega_t a k =
  k.wt *. a /. (a +. k.theta_nr)

let omega_q a k q =
  k.wq *. a /. ((a +. k.theta_nr) *. (1. +. (q /. k.kKq) ** k.hq))

let nu_t c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a 1.0) in
  rate /. k.nt

let nu_m c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a 1.0) in
  rate /. k.nm

let nu_q c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a 1.0) in
  rate /. k.nq

let nu_r c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a 1.0) in
  rate /. k.nr

(* Kept for reference purpose *)
(* let compute_mM' r t m q u c_r c_t c_m c_q c_u k = *)
(*   k.nr *. r +. k.nt *. t +. k.nm *. m +. k.nq *. q +. k.nu *. u +. k.nr *. (c_r +. c_t +. c_m +. c_q +. c_u) *)

(* let compute_lambda a r e_t e_m q c_r c_t c_m c_q k usemMparam = *)
(*   let mM = *)
(*     if usemMparam then *)
(*       k.mM *)
(*     else *)
(*       compute_mM' r e_t e_m q c_r c_t c_m c_q k *)
(*   in *)
(*   let vRt   = c_r +. c_t +. c_m +. c_q in *)
(*   let gamma = hill k.gamma_max k.kK_gamma a 1.0 in *)
(*   gamma *. vRt /. mM *)

(* Compute protein mass *)
let compute_mM
    { r   = { prot = r;   bound_mrna = c_r };
      e_t = { prot = e_t; bound_mrna = c_t };
      e_m = { prot = e_m; bound_mrna = c_m };
      q   = { prot = q;   bound_mrna = c_q }
    }
    k =
  k.nr *. r +.  k.nt *. e_t +.  k.nm *. e_m +.  k.nq *. q +.  k.nr *. (c_r +. c_t +. c_m +. c_q)

(* the useM parameter is an artefact of matlabish code - to be removed. *)
let compute_lambda { a; r; e_t; e_m; q } k useM =
  let mM = k.mM
    (* if useM then *)
    (*   k.mM *)
    (* else *)
    (*   compute_mM cell k *)
  in
  let vRt   = r.bound_mrna +. e_t.bound_mrna +. e_m.bound_mrna +. q.bound_mrna in
  let gamma = hill k.gamma_max k.kK_gamma a 1.0 in
  gamma *. vRt /. mM

(* Conversion from the CVODE state array to a nice record.   *)
let cell_of_state (state : Sundials.RealArray.t) k =
  let s_i       = state.{0} in
  let a         = state.{1} in
  let r         = state.{2} in
  let e_t       = state.{3} in
  let e_m       = state.{4} in
  let q         = state.{5} in
  let m_r       = state.{6} in
  let m_t       = state.{7} in
  let m_m       = state.{8} in
  let m_q       = state.{9} in
  let c_r       = state.{10} in
  let c_t       = state.{11} in
  let c_m       = state.{12} in
  let c_q       = state.{13} in

  {
    s_i;
    a;
    r = {
      prot       = r;
      free_mrna  = m_r;
      bound_mrna = c_r;

      nu         = (nu_r c_r a k);
      omega      = (omega_r a k);
    };
    e_t = {
      prot       = e_t;
      free_mrna  = m_t;
      bound_mrna = c_t;

      nu         = (nu_t c_t a k);
      omega      = (omega_t a k);
    };
    e_m = {
      prot       = e_m;
      free_mrna  = m_m;
      bound_mrna = c_m;

      nu         = (nu_m c_m a k);
      omega      = (omega_m a k);
    };
    q = {
      prot       = q;
      free_mrna  = m_q;
      bound_mrna = c_q;

      nu         = (nu_q c_q a k);
      omega      = (omega_q a k q);
    }
  }

let ddt_prot { prot; nu } lambda =
  nu -. lambda *. prot

(* ribosome protein production *)    
let ddt_prot_r r e_t e_m q lambda k =
      r.nu   -. lambda *. r.prot
  +. (r.nu   -. k.k_b *. r.prot *. r.free_mrna   +. k.k_u *. r.bound_mrna)
  +. (e_t.nu -. k.k_b *. r.prot *. e_t.free_mrna +. k.k_u *. e_t.bound_mrna)
  +. (e_m.nu -. k.k_b *. r.prot *. e_m.free_mrna +. k.k_u *. e_m.bound_mrna)
  +. (q.nu   -. k.k_b *. r.prot *. q.free_mrna   +. k.k_u *. q.bound_mrna)

let ddt_free_mrna ribosome { prot; free_mrna; bound_mrna; omega; nu } lambda k =
  omega -. (lambda +. k.dm) *. free_mrna +. nu -. k.k_b *. ribosome.prot *. free_mrna +. k.k_u *. bound_mrna

let ddt_bound_mrna ribosome { bound_mrna; free_mrna; nu } lambda k =
  -. lambda *. bound_mrna
  +. k.k_b  *. ribosome.prot *. free_mrna
  -. k.k_u *. bound_mrna
  -. nu
  


(********************************************************)
(* Tests implementation - UGLY *)
(********************************************************)    

let epsilon = 1e-3

let nearzero x scale = (abs_float x) < (epsilon *. scale)

(* let test1 cell lambda k = *)
(*   let x = *)
(*     lambda *. k.mM /. k.ns *)
(*   in *)
(*   let y = *)
(*     cell.e_m.prot *. k.vm /. (1. +. k.kKm /. cell.s_i) *)
(*   in *)
(*   let z = *)
(*     cell.e_t.prot *. k.vt /. (1. +. k.kKt /. k.s) *)
(*   in *)
(*   let _ = printf "%f %f %f\n%!" x y z in *)
(*   let scale =  (abs_float x) +. (abs_float y) +. (abs_float z) in *)
(*   (nearzero (x -. y) scale) && *)
(*   (nearzero (x -. z) scale) && *)
(*   (nearzero (y -. z) scale) *)

(* let test2 cell k lambda = *)
(*   let gamma = hill k.gamma_max k.kK_gamma cell.a 1.0 in *)
(*   let ratio = *)
(*     1. /. (cell.r.bound_mrna +. *)
(*            cell.e_t.bound_mrna +. *)
(*            cell.e_m.bound_mrna +. *)
(*            cell.q.bound_mrna +. *)
(*            cell.r.prot) *)
(*   in *)
(*   let (betar, betat, betam, betaq) = (cell.r.bound_mrna   *. ratio, *)
(*                                       cell.e_t.bound_mrna *. ratio, *)
(*                                       cell.e_m.bound_mrna *. ratio, *)
(*                                       cell.q.bound_mrna   *. ratio) in *)
(*   let ur  = lambda *. k.mM *. betar in *)
(*   let ut  = lambda *. k.mM *. betat in *)
(*   let um  = lambda *. k.mM *. betam in *)
(*   let uq  = lambda *. k.mM *. betaq in *)

(*   let vr = *)
(*     let a  = gamma *. k.wr /. (1. +. k.theta_r /. cell.a) *)
(*     and b =  *)
(*       let b1 = lambda *)
(*       and b2 = lambda +. k.dm *)
(*       and b3 = *)
(*         let num   = lambda +. k.k_u +. gamma /. k.nr *)
(*         and denum = k.k_b *. cell.r.prot *)
(*         in *)
(*         num /. denum *)
(*       in *)
(*       b1 +. b2 *. b3 *)
(*     in *)
(*     a /. b *)
(*   in *)
(*   let vt = *)
(*     let a = gamma *. k.wt /. (1. +. k.theta_nr /. cell.a)  *)
(*     and b = *)
(*       let b1 = lambda *)
(*       and b2 = lambda +. k.dm *)
(*       and b3 = *)
(*         let num   = lambda +. k.k_u +. gamma /. k.nt *)
(*         and denum = k.k_b *. cell.e_t.prot *)
(*         in *)
(*         num /. denum *)
(*       in *)
(*       b1 +. b2 *. b3 *)
(*     in *)
(*     a /. b *)

(*   in *)
(*   let _ = printf "%f %f %f %f\n%!" ur vr ut vt in *)
(*   let scale =  (abs_float ur) in *)
(*   let scale' =  (abs_float ut) in *)
(*   (nearzero (ur -. vr) scale) && *)
(*   (nearzero (ut -. vt) scale') *)


  

(********************************************************)
(* Step function *)
(********************************************************)

let step k (time : float) (x : Sundials.RealArray.t) (dxdt : Sundials.RealArray.t) =

  let cell   = cell_of_state x k in

  let lambda = compute_lambda cell k true in

  let vimp = cell.e_t.prot *. (hill k.vt k.kKt k.s 1.0) in
  let vcat = cell.e_m.prot *. (hill k.vm k.kKm cell.s_i 1.0) in

  let ddts_i     = vimp -. vcat -. lambda *. cell.s_i in

  let ddta       =   k.ns *. vcat
                   -. (k.nr *. cell.r.nu)
                   -. (k.nt *. cell.e_t.nu)
                   -. (k.nm *. cell.e_m.nu)
                   -. (k.nq *. cell.q.nu)
                   -. lambda *. cell.a in

  let ddtr   = ddt_prot_r cell.r cell.e_t cell.e_m cell.q lambda k in
  let ddte_t = ddt_prot cell.e_t lambda in
  let ddte_m = ddt_prot cell.e_m lambda in
  let ddtq   = ddt_prot cell.q lambda in

  let ddtm_r = ddt_free_mrna cell.r cell.r lambda k in
  let ddtm_t = ddt_free_mrna cell.r cell.e_t lambda k in
  let ddtm_m = ddt_free_mrna cell.r cell.e_m lambda k in
  let ddtm_q = ddt_free_mrna cell.r cell.q lambda k in

  let ddtc_r = ddt_bound_mrna cell.r cell.r lambda k in
  let ddtc_t = ddt_bound_mrna cell.r cell.e_t lambda k in
  let ddtc_m = ddt_bound_mrna cell.r cell.e_m lambda k in
  let ddtc_q = ddt_bound_mrna cell.r cell.q lambda k in

  (* if time > 2800.0 then *)
  (*   assert (test2 cell k lambda) *)
  (* else (); *)

  dxdt.{0} <- ddts_i; 
  dxdt.{1} <- ddta; 
  dxdt.{2} <- ddtr; 
  dxdt.{3} <- ddte_t; 
  dxdt.{4} <- ddte_m; 
  dxdt.{5} <- ddtq; 
  dxdt.{6} <- ddtm_r; 
  dxdt.{7} <- ddtm_t; 
  dxdt.{8} <- ddtm_m; 
  dxdt.{9} <- ddtm_q; 
  dxdt.{10} <- ddtc_r; 
  dxdt.{11} <- ddtc_t; 
  dxdt.{12} <- ddtc_m; 
  dxdt.{13} <- ddtc_q

(********************************************************)
(* Initialisation *)
(********************************************************)

  
let initial_parameters = {
  s        = 1e4; (* nutrient concentration *)
  dm       = 0.1; (* mRNA degradation rate *)
  ns       = 0.5;
  nr       = 7459.0;
  nt       = 300.0;
  nm       = 300.0;
  nq       = 300.0;
  gamma_max = 1260.0;
  kK_gamma  = 7.0;
  vt        = 726.0;
  kKt       = 1000.0;
  vm        = 5800.0;
  kKm       = 1000.0;
  wr        = 930.0;
  wt        = 4.14;
  wm        = 4.14;
  wq        = 948.93;
  theta_r   = 426.87;
  theta_nr  = 4.38;
  kKq       = 152219.0;
  hq        = 4.0;
  k_b       = 1.0;
  k_u       = 1.0;
  mM        = 1e8;
}

(* TODO: give proper names to the contents of the cells of this array. *)

let initial_state =
  [| 0.0; 0.0; 1.0; 1.0; 1.0; 0.0; 0.0;
     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;
  |]

let result_state =
  [| 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;
     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;
  |]

let stop_time = 3000.0

(* Integrator setup *)
let x0 = Nvector_serial.wrap (Sundials.RealArray.of_array initial_state)

let xt = Nvector_serial.wrap (Sundials.RealArray.of_array result_state)

let tolerances = Cvode.SStolerances(1e-3, 1e-6)

(* Configure the solver for a stiff (Cvode.BDF), dense problem. See Cvode.mli *)
let s = Cvode.init Cvode.BDF (Cvode.Newton (Cvode.Dls.dense ())) tolerances (step initial_parameters)  0. x0

let _ = Cvode.set_stop_time s stop_time

(* let compute_proteic_mass x k = *)
(*   compute_mM' x.{2} x.{3} x.{4} x.{5} x.{10} x.{11} x.{12} x.{13} k *)

let fd, print_state = fprint_vec "state"
let print_mass x  = fprint "mass" x

(* Integration *)
let _ =
  for i = 1 to (int_of_float stop_time) do
    let (_, _) = Cvode.solve_normal s (float i) xt in
    print_state (Nvector.unwrap xt)
  done;
  close_out fd


