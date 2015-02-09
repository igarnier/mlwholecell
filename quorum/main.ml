(* Interconnecting cells
   sensing AHL in neighbouring cells
   if hill(AHL) is present, this triggers the promoter q
   when hill(promoter q) is present,
   q gets produced *)

module RealArray = Sundials.RealArray
(* module Roots = Sundials.Roots *)

type params = {
  s         : float; (* nutrient concentration *)
  dm        : float; (* mRNA degradation rate *)
  ns        : float; (* nutrient efficiency during transcription *)
  nr        : float; (* protein lengths *)
  nt        : float;
  nm        : float;
  nq        : float;
  nAHL      : float;
  nLuxL     : float;
  nLuxR     : float;
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
  wAHL      : float;
  wLuxL     : float;
  wLuxR     : float;
  theta_r   : float; (* ribosome transcription threshold  *)
  theta_nr  : float;
  kKq       : float; (* q-qutoinhibition threshold *)
  hq        : float; (* q-qutoinhibition Hill coef. *)
  k_b       : float; (* mRNA-ribosome binding rate *)
  k_u       : float; (* mRNA-ribosome unbinding rate *)
  k_b_lux_ahl : float; (* LuxR-AHL binding rate *)
  k_u_lux_ahl : float; (* LuxR-AHL unbinding rate *)
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
  s_i     : float;
  a       : float;
  r       : protein;
  e_t     : protein;
  e_m     : protein;
  q       : protein;
  ahl     : protein;
  luxl    : protein;
  luxr    : protein;
  lux_ahl : float    (* luxr-AHL complex *)
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

let fprint_vec pattern file =
  let fd = open_out file in
  fun (x : Sundials.RealArray.t option) ->
    match x with
    | None -> close_out fd
    | Some data ->
      (* let bigarr = Nvector.unwrap data in *)
      for i = 0 to Array.length pattern - 1 do
        fprintf fd "%f " data.{pattern.(i)}
      done;
      fprintf fd "\n"

(********************************************************)
(* Helper functions *)
(********************************************************)

(* Variables naming *)
let var_s_i    = 0
let var_a      = 1
let var_r      = 2
let var_e_t    = 3
let var_e_m    = 4
let var_q      = 5
let var_m_r    = 6
let var_m_t    = 7
let var_m_m    = 8
let var_m_q    = 9
let var_c_r    = 10
let var_c_t    = 11
let var_c_m    = 12
let var_c_q    = 13
let var_m_ahl  = 14
let var_c_ahl  = 15
let var_luxl   = 16
let var_m_luxl = 17
let var_c_luxl = 18
let var_luxr   = 19
let var_m_luxr = 20
let var_c_luxr = 21
let var_luxr_ahl = 22

(* Hill function *)
let hill vmax k_half_rate x =
  vmax *. x /. (k_half_rate +. x)

(* Transcription rates *)

let omega_r a k =
  k.wr *. a /. (a +. k.theta_r)

let omega_m a k =
  k.wm *. a /. (a +. k.theta_nr)

let omega_t a k =
  k.wt *. a /. (a +. k.theta_nr)

let omega_q a k q =
  k.wq *. a /. ((a +. k.theta_nr) *. (1. +. (q /. k.kKq) ** k.hq))

let omega_AHL a luxl k =
  k.wAHL *. a /. (a +. k.theta_nr) *. luxl /. (luxl +. k.theta_nr)

let omega_LuxL a lux_ahl k =
  k.wLuxL *. a /. (a +. k.theta_nr) *. lux_ahl /. (lux_ahl +. k.theta_nr *. 100.0)

let omega_LuxR a k =
  k.wLuxR *. a /. (a +. k.theta_nr)

(* Translation rates *)

let nu_t c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nt

let nu_m c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nm

let nu_q c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nq

let nu_r c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nr

let nu_AHL c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nAHL

let nu_LuxL c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nLuxL

let nu_LuxR c_x a k =
  let rate = c_x *. (hill k.gamma_max k.kK_gamma a) in
  rate /. k.nLuxR

let compute_mM' r t m q ahl luxl luxr c_r c_t c_m c_q c_AHL c_LuxL c_LuxR lux_ahl k =
  k.nr *. r 
  +. k.nt *. t
  +. k.nm *. m 
  +. k.nq *. q 
  +. k.nAHL *. ahl 
  +. k.nLuxL *. luxl
  +. k.nLuxR *. luxr
  +. (k.nLuxR +. k.nAHL) *. lux_ahl
  +. k.nr *. (c_r +. c_t +. c_m +. c_q +. c_AHL +. c_LuxL +. c_LuxR)

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

let compute_mM
    { r    = { prot = r;    bound_mrna = c_r; _ };
      e_t  = { prot = e_t;  bound_mrna = c_t; _ };
      e_m  = { prot = e_m;  bound_mrna = c_m; _ };
      q    = { prot = q;    bound_mrna = c_q; _ }; 
      ahl  = { prot = ahl;  bound_mrna = c_ahl; _ };
      luxl = { prot = luxl; bound_mrna = c_luxl; _ };
      luxr = { prot = luxr; bound_mrna = c_luxr; _ };
      lux_ahl;
      _ } k =
  k.nr *. r  
  +. k.nt *. e_t
  +. k.nm *. e_m 
  +. k.nq *. q 
  +. k.nAHL *. ahl 
  +. k.nLuxL *. luxl
  +. k.nLuxR *. luxr
  +. (k.nLuxR +. k.nAHL) *. lux_ahl
  +. k.nr *. (c_r +. c_t +. c_m +. c_q +. c_ahl +. c_luxl +. c_luxr)

let compute_lambda ({ a; r; e_t; e_m; q; ahl; luxl; luxr; _ } as cell) k useM =
  let mM =
    if useM then
      k.mM
    else
      compute_mM cell k
  in
  let vRt   = r.bound_mrna +. e_t.bound_mrna +. e_m.bound_mrna +. q.bound_mrna +. ahl.bound_mrna +. luxl.bound_mrna +. luxr.bound_mrna in
  let gamma = hill k.gamma_max k.kK_gamma a in
  gamma *. vRt /. mM

let cell_of_state base_index ahl_input (state : Sundials.RealArray.t) k =
  let s_i       = state.{base_index + 0} in
  let a         = state.{base_index + 1} in
  let r         = state.{base_index + 2} in
  let e_t       = state.{base_index + 3} in
  let e_m       = state.{base_index + 4} in
  let q         = state.{base_index + 5} in
  let m_r       = state.{base_index + 6} in
  let m_t       = state.{base_index + 7} in
  let m_m       = state.{base_index + 8} in
  let m_q       = state.{base_index + 9} in
  let c_r       = state.{base_index + 10} in
  let c_t       = state.{base_index + 11} in
  let c_m       = state.{base_index + 12} in
  let c_q       = state.{base_index + 13} in

  let ahl       = state.{ahl_input} in
  let m_ahl     = state.{base_index + 14} in
  let c_ahl     = state.{base_index + 15} in

  let luxl      = state.{base_index + 16} in
  let m_luxl    = state.{base_index + 17} in
  let c_luxl    = state.{base_index + 18} in

  let luxr      = state.{base_index + 19} in
  let m_luxr    = state.{base_index + 20} in
  let c_luxr    = state.{base_index + 21} in

  let lux_ahl   = state.{base_index + 22} in

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
    };
    ahl = {
      prot       = ahl;
      free_mrna  = m_ahl;
      bound_mrna = c_ahl;
      
      nu         = (nu_AHL c_ahl a k);
      omega      = (omega_AHL a luxl k);      
    };
    luxl = {
      prot       = luxl;
      free_mrna  = m_luxl;
      bound_mrna = c_luxl;
      
      nu         = (nu_LuxL c_luxl a k);
      omega      = (omega_LuxL a lux_ahl k);
    };
    luxr = {
      prot       = luxr;
      free_mrna  = m_luxr;
      bound_mrna = c_luxr;
      
      nu         = (nu_LuxR c_luxr a k);
      omega      = (omega_LuxR a k);      
    };
    lux_ahl
  }

let ddt_prot { prot; nu; _ } lambda =
  nu -. lambda *. prot
    
let ddt_prot_r r e_t e_m q ahl luxl luxr lambda k =
      r.nu   -. lambda *. r.prot
  +. (r.nu   -. k.k_b *. r.prot *. r.free_mrna   +. k.k_u *. r.bound_mrna)
  +. (e_t.nu -. k.k_b *. r.prot *. e_t.free_mrna +. k.k_u *. e_t.bound_mrna)
  +. (e_m.nu -. k.k_b *. r.prot *. e_m.free_mrna +. k.k_u *. e_m.bound_mrna)
  +. (q.nu   -. k.k_b *. r.prot *. q.free_mrna   +. k.k_u *. q.bound_mrna)
  +. (ahl.nu   -. k.k_b *. r.prot *. ahl.free_mrna +. k.k_u *. ahl.bound_mrna)
  +. (luxl.nu   -. k.k_b *. r.prot *. luxl.free_mrna +. k.k_u *. luxl.bound_mrna)
  +. (luxr.nu   -. k.k_b *. r.prot *. luxr.free_mrna +. k.k_u *. luxr.bound_mrna)

let ddt_free_mrna ribosome { free_mrna; bound_mrna; omega; nu; _ } lambda k =
  omega                                   (* Transcription rate *)
  +. nu                                   (* After transcription is done, the bound mrna becomes free again *)
  +. k.k_u *. bound_mrna                  (* Spontaneous unbinding *)
  -. (lambda +. k.dm) *. free_mrna
  -. k.k_b *. ribosome.prot *. free_mrna 

let ddt_bound_mrna ribosome { bound_mrna; free_mrna; nu; _ } lambda k =
  -. lambda *. bound_mrna
  +. k.k_b  *. ribosome.prot *. free_mrna
  -. k.k_u *. bound_mrna
  -. nu

let ddt_luxr_ahl ahl luxr lux_ahl lambda k =
  -. lambda *. lux_ahl                            (* Dilution term *)
  +. ahl.prot *. luxr.prot *. k.k_b_lux_ahl       (* binding ahl and luxr *)
  -. lux_ahl *. k.k_u_lux_ahl                     (* unbinding ahl and luxr *)

let ddt_prot_ahl { prot; nu; _ } luxr lux_ahl lambda k =
  nu                                              (* outcome from translation *)
  -. lambda *. prot                               (* dilution *)
  -. prot *. luxr.prot *. k.k_b_lux_ahl           (* binding ahl and luxr *)
  +. lux_ahl *. k.k_u_lux_ahl                     (* unbinding ahl and luxr *)

let ddt_prot_luxr { prot; nu; _ } ahl lux_ahl lambda k =
  nu                                              (* outcome from translation *)
  -. lambda *. prot                               (* dilution *)
  -. prot *. ahl.prot *. k.k_b_lux_ahl            (* binding ahl and luxr *)
  +. lux_ahl *. k.k_u_lux_ahl                     (* unbinding ahl and luxr *)


(********************************************************)
(* Tests implementation *)
(********************************************************)    

let epsilon = 1e-3

let nearzero x scale = (abs_float x) < (epsilon *. scale)

(*
let test1 cell lambda k =
  let x =
    lambda *. k.mM /. k.ns
  in
  let y =
    cell.e_m.prot *. k.vm /. (1. +. k.kKm /. cell.s_i)
  in
  let z =
    cell.e_t.prot *. k.vt /. (1. +. k.kKt /. k.s)
  in
  let _ = printf "%f %f %f\n%!" x y z in
  let scale =  (abs_float x) +. (abs_float y) +. (abs_float z) in
  (nearzero (x -. y) scale) &&
  (nearzero (x -. z) scale) &&
  (nearzero (y -. z) scale)

let test2 cell k lambda =
  let gamma = hill k.gamma_max k.kK_gamma cell.a 1.0 in
  let ratio =
    1. /. (cell.r.bound_mrna +.
           cell.e_t.bound_mrna +.
           cell.e_m.bound_mrna +.
           cell.q.bound_mrna +.
           cell.r.prot)
  in
  let (betar, betat, betam, betaq) = (cell.r.bound_mrna   *. ratio,
                                      cell.e_t.bound_mrna *. ratio,
                                      cell.e_m.bound_mrna *. ratio,
                                      cell.q.bound_mrna   *. ratio) in
  let ur  = lambda *. k.mM *. betar in
  let ut  = lambda *. k.mM *. betat in
  let um  = lambda *. k.mM *. betar in
  let uq  = lambda *. k.mM *. betar in

  let vr =
    let a  = gamma *. k.wr /. (1. +. k.theta_r /. cell.a)
    and b = 
      let b1 = lambda
      and b2 = lambda +. k.dm
      and b3 =
        let num   = lambda +. k.k_u +. gamma /. k.nr
        and denum = k.k_b *. cell.r.prot
        in
        num /. denum
      in
      b1 +. b2 *. b3
    in
    a /. b
  in
  let vt =
    let a = gamma *. k.wt /. (1. +. k.theta_nr /. cell.a) 
    and b =
      let b1 = lambda
      and b2 = lambda +. k.dm
      and b3 =
        let num   = lambda +. k.k_u +. gamma /. k.nt
        and denum = k.k_b *. cell.e_t.prot
        in
        num /. denum
      in
      b1 +. b2 *. b3
    in
    a /. b

  in
  let _ = printf "%f %f %f %f\n%!" ur vr ut vt in
  let scale =  (abs_float ur) in
  let scale' =  (abs_float ut) in
  (nearzero (ur -. vr) scale) &&
  (nearzero (ut -. vt) scale')
*)

(********************************************************)
(* Diffusion-related stuff *)
(********************************************************)

let index_of_coords x y ysize =
  x * ysize + y

let coords_of_index index ysize =
  (index / ysize, index mod ysize)

let indexing xsize ysize =
  let a = Array.make_matrix xsize ysize 0 in
  for x = 0 to xsize - 1 do
    for y = 0 to ysize - 1 do
      a.(x).(y) <- index_of_coords x y ysize
    done
  done;
  a
      
let grad (state : Sundials.RealArray.t) (index : int array array) x y =
  let w = state.{index.(x - 1).(y)}
  and e = state.{index.(x + 1).(y)}
  and s = state.{index.(x).(y - 1)}
  and n = state.{index.(x).(y + 1)}
  in 
  w +. e +. s +. n -. 4. *. state.{index.(x).(y)}      


(********************************************************)
(* Memory layout of the state *)
(********************************************************)

(* The first tissue_X * tissue_Y indices in the state
   are dedicated to storing the AHL concentrations. After
   that, the variables for the cells which are present
   (according to matrix [presence]) are stored contiguously.
   The value of the [presence] matrix at (x, y) is either:
   . -1 if no cells are present at (x, y)
   . some offset index p in the state, corresponding to 
     where the cell is stored.
   The [gradient] auxilliary array is dedicated to the
   storage of the gradient, recomputed at each step. This
   avoids reallocation. [diffusion] is the diffusion 
   coefficient. It should be fitted relativelw to a
   particular unit grid length.
*)

(********************************************************)
(* Step functions *)
(********************************************************)

let step_cell output_index base_index cell k (dxdt : Sundials.RealArray.t) =
  let lambda = compute_lambda cell k true in

  (* let _      = printf "lambda %f\n" lambda in *)

  let vimp = cell.e_t.prot *. (hill k.vt k.kKt k.s) in
  let vcat = cell.e_m.prot *. (hill k.vm k.kKm cell.s_i) in

  let ddts_i     = vimp -. vcat -. lambda *. cell.s_i in

  let ddta       =   k.ns *. vcat
                   -. (k.nr *. cell.r.nu)
                   -. (k.nt *. cell.e_t.nu)
                   -. (k.nm *. cell.e_m.nu)
                   -. (k.nq *. cell.q.nu)
                   -. (k.nAHL *. cell.ahl.nu)
                   -. (k.nLuxL *. cell.luxl.nu)
                   -. (k.nLuxR *. cell.luxr.nu)
                   -. lambda *. cell.a in

  let ddtr     = ddt_prot_r cell.r cell.e_t cell.e_m cell.q cell.ahl cell.luxl cell.luxr lambda k in
  let ddte_t   = ddt_prot cell.e_t lambda in
  let ddte_m   = ddt_prot cell.e_m lambda in
  let ddtq     = ddt_prot cell.q lambda in
  let ddtahl   = ddt_prot_ahl cell.ahl cell.luxr cell.lux_ahl lambda k in
  let ddtluxl  = ddt_prot cell.luxl lambda in
  let ddtluxr  = ddt_prot_luxr cell.luxr cell.ahl cell.lux_ahl lambda k in

  let ddtm_r    = ddt_free_mrna cell.r cell.r lambda k in
  let ddtm_t    = ddt_free_mrna cell.r cell.e_t lambda k in
  let ddtm_m    = ddt_free_mrna cell.r cell.e_m lambda k in
  let ddtm_q    = ddt_free_mrna cell.r cell.q lambda k in
  let ddtm_ahl  = ddt_free_mrna cell.r cell.ahl lambda k in
  let ddtm_luxl = ddt_free_mrna cell.r cell.luxl lambda k in
  let ddtm_luxr = ddt_free_mrna cell.r cell.luxr lambda k in

  let ddtc_r    = ddt_bound_mrna cell.r cell.r lambda k in
  let ddtc_t    = ddt_bound_mrna cell.r cell.e_t lambda k in
  let ddtc_m    = ddt_bound_mrna cell.r cell.e_m lambda k in
  let ddtc_q    = ddt_bound_mrna cell.r cell.q lambda k in  
  let ddtc_ahl  = ddt_bound_mrna cell.r cell.ahl lambda k in
  let ddtc_luxl = ddt_bound_mrna cell.r cell.luxl lambda k in
  let ddtc_luxr = ddt_bound_mrna cell.r cell.luxr lambda k in

  let ddt_luxr_ahl = ddt_luxr_ahl cell.ahl cell.luxr cell.lux_ahl lambda k in

  dxdt.{base_index + 0} <- ddts_i; 
  dxdt.{base_index + 1} <- ddta; 
  dxdt.{base_index + 2} <- ddtr; 
  dxdt.{base_index + 3} <- ddte_t; 
  dxdt.{base_index + 4} <- ddte_m; 
  dxdt.{base_index + 5} <- ddtq; 
  dxdt.{base_index + 6} <- ddtm_r; 
  dxdt.{base_index + 7} <- ddtm_t; 
  dxdt.{base_index + 8} <- ddtm_m; 
  dxdt.{base_index + 9} <- ddtm_q; 
  dxdt.{base_index + 10} <- ddtc_r; 
  dxdt.{base_index + 11} <- ddtc_t; 
  dxdt.{base_index + 12} <- ddtc_m; 
  dxdt.{base_index + 13} <- ddtc_q;
  dxdt.{output_index}    <- ddtahl;
  dxdt.{base_index + 14} <- ddtm_ahl;
  dxdt.{base_index + 15} <- ddtc_ahl;
  
  dxdt.{base_index + 16} <- ddtluxl;
  dxdt.{base_index + 17} <- ddtm_luxl;
  dxdt.{base_index + 18} <- ddtc_luxl;
  
  dxdt.{base_index + 19} <- ddtluxr;
  dxdt.{base_index + 20} <- ddtm_luxr;
  dxdt.{base_index + 21} <- ddtc_luxr;
  dxdt.{base_index + 22} <- ddt_luxr_ahl
  


(* First, compute cells evolution. Then, update diffusion state with outcome. *)  
let step
    (xsize : int)
    (ysize : int)
    (diffusion : float)
    (index : int array array)
    (gradient : Sundials.RealArray.t)
    (presence : int array array)
    (k : params)
    (time : float)
    (state : Sundials.RealArray.t)
    (dxdt : Sundials.RealArray.t) =
  (* Warning: will segfault if cells are present on the boundary. *)
  (* Reset dxdt *)
  for i = 0 to Sundials.RealArray.length dxdt - 1 do
    dxdt.{i} <- 0.0
  done;
  (* Implement absorbing boundary *)
  for x = 0 to xsize - 1 do
    state.{index.(x).(0)} <- 0.0;
    state.{index.(x).(ysize - 1)} <- 0.0;
  done;
  for y = 0 to ysize - 1 do
    state.{index.(0).(y)} <- 0.0;
    state.{index.(xsize - 1).(y)} <- 0.0
  done;
  for x = 0 to xsize - 1 do
    for y = 0 to ysize - 1 do
      let cell_idx = presence.(x).(y) in
      if cell_idx >= 0 then
        begin
          let cell = cell_of_state cell_idx index.(x).(y) state k in
          step_cell index.(x).(y) cell_idx cell k dxdt
          (* Place produced AHL outside of the cell *)
          (* let ahl_produced = dxdt.{cell_idx + var_ahl} *. 0.25 in *)
          (* dxdt.{cell_idx + var_ahl} <- 0.0; *)
          (* dxdt.{index.(x - 1).(y)} <- dxdt.{index.(x - 1).(y)} +. ahl_produced; *)
          (* dxdt.{index.(x + 1).(y)} <- dxdt.{index.(x + 1).(y)} +. ahl_produced; *)
          (* dxdt.{index.(x).(y - 1)} <- dxdt.{index.(x).(y - 1)} +. ahl_produced; *)
          (* dxdt.{index.(x).(y + 1)} <- dxdt.{index.(x).(y + 1)} +. ahl_produced *)
        end
      else ()
    done
  done;
  (* Update gradient *)
  for x = 1 to xsize - 2 do
    for y = 1 to ysize - 2 do
      gradient.{index.(x).(y)} <- (grad state index x y)
    done
  done;
  (* Update dxdt *)
  for x = 1 to xsize - 2 do
    for y = 1 to ysize - 2 do
      let delta = grad gradient index x y in
      dxdt.{index.(x).(y)} <- dxdt.{index.(x).(y)}  -. diffusion *. delta
    done
  done


(********************************************************)
(* Graphical display *)
(********************************************************)

let interpolate (r1,g1,b1) (r2,g2,b2) f =
  let deltar = float (r2 - r1) in
  let deltag = float (g2 - g1) in
  let deltab = float (b2 - b1) in
  let r     = int_of_float (f *. deltar +. (float r1))
  and g     = int_of_float (f *. deltag +. (float g1))
  and b     = int_of_float (f *. deltab +. (float b1)) in
  Graphics.rgb r g b

(* Note that the palette is skewed towards low values. *)
let palette =
  let palette = Array.make 256 Graphics.black in
  let x1      = 85.0 in
  let x2      = 85.0 *. 2.0 in
  for i = 0 to (int_of_float x1) - 1 do
    let ratio = (float i) /. x1 in
    palette.(i) <- interpolate (0, 0, 0) (255, 0, 0) ratio
  done;
  for i = (int_of_float x1) to (int_of_float x2) - 1 do
    let ratio = (float i -. x1) /. (x2 -. x1) in
    palette.(i) <- interpolate (255, 0, 0) (0, 255, 0) ratio
  done;
  for i = (int_of_float x2) to 255 do
    let ratio = (float i -. x2) /. (256.0 -. x2) in
    palette.(i) <- interpolate (0, 255, 255) (0, 0, 255) ratio
  done;
  palette

(* make this a parameter of the program  *)
let max_value = float_of_string Sys.argv.(2)

let redresse x = -. 1. /. (30. *. x +. 1.) +. 1.

let float_to_color f =
  let f = min max_value f in
  let r = redresse (f /. max_value) in
  palette.(int_of_float (r *. 255.))

let display_cells tissue_X tissue_Y square_size presence tissue var =
  Graphics.set_color Graphics.black;
  Graphics.fill_rect
    0 0
    (tissue_X * square_size)
    (tissue_X * square_size);
  for x = 0 to tissue_X - 1 do
    for y = 0 to tissue_Y - 1 do
      let i = presence.(x).(y) in
      if i >= 0 then
        let cell = abs_float tissue.{i + var} in
        let clr  = float_to_color cell in
        Graphics.set_color clr;
        Graphics.fill_rect
          (x * square_size) 
          (y * square_size)
          square_size
          square_size
      else ()
    done
  done;
  Graphics.synchronize ()

let display_diffusion tissue_X tissue_Y square_size index state =
  Graphics.set_color Graphics.black;
  Graphics.fill_rect
    0 0
    (tissue_X * square_size)
    (tissue_X * square_size);
  for x = 0 to tissue_X - 1 do
    for y = 0 to tissue_Y - 1 do
      let cell = abs_float state.{index.(x).(y)} in
      let clr  = float_to_color cell in
      Graphics.set_color clr;
      Graphics.fill_rect
        (x * square_size) 
        (y * square_size)
        square_size
        square_size
    done
  done


let display_info time =
  Graphics.set_color Graphics.white;
  Graphics.moveto 30 30;
  Graphics.draw_string (string_of_int time)


(********************************************************)
(* Initialisation *)
(********************************************************)

(* Space size *)
let tissue_X = 50
let tissue_Y = 50

(* Number of array slots dedicated to diffusion *)
let diffusion_size  = tissue_X * tissue_Y

let initial_AHL_level = 0.001
  
let initial_parameters = {
  s         = 1e4; (* nutrient concentration *)
  dm        = 0.1; (* mRNA degradation rate *)
  ns        = 0.5;
  nr        = 7459.0;
  nt        = 300.0;
  nm        = 300.0;
  nq        = 300.0;
  nAHL      = 300.0;
  nLuxL     = 300.0;
  nLuxR     = 300.0;
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
  wAHL      = 0.1; (* TODO, max transcription rate for AHL *)
  wLuxL     = 0.1; (* TODO, max transcription rate for LuxL *)
  wLuxR     = 0.1; (* TODO, max transcription rate for LuxR *)
  theta_r   = 426.87;
  theta_nr  = 4.38;
  kKq       = 152219.0;
  hq        = 4.0;
  k_b       = 1.0;
  k_u       = 1.0;
  k_b_lux_ahl = 1.0;
  k_u_lux_ahl = 1.0;
  mM        = 1e8;
}


(* USER-DEFINED initial state for all cells. *)
let initial_state_per_cell =
  let a = Array.make 23 0.0 in
  a.(var_s_i) <- 0.0;
  a.(var_a)   <- 0.0;
  a.(var_r)   <- 1.0;
  a.(var_e_t) <- 1.0;
  a.(var_e_m) <- 1.0;
  a.(var_q)   <- 0.0;
  a.(var_m_r) <- 0.0;
  a.(var_m_t) <- 0.0;
  a.(var_m_m) <- 0.0;
  a.(var_m_q) <- 0.0;
  a.(var_c_r) <- 0.0;
  a.(var_c_t) <- 0.0;
  a.(var_c_m) <- 0.0;
  a.(var_c_q) <- 0.0;
  a.(var_m_ahl) <- 0.0;
  a.(var_c_ahl) <- 0.0;
  a.(var_luxl) <- 0.0;
  a.(var_m_luxl) <- 0.0;
  a.(var_c_luxl) <- 0.0;
  a.(var_luxr) <- 0.0;
  a.(var_m_luxr) <- 0.0;
  a.(var_c_luxr) <- 0.0;
  a.(var_luxr_ahl) <- 0.0;
  a

let variables_count = Array.length initial_state_per_cell

(* USER-DEFINED boolean presence matrix *)
(* Initialise with just one cell in the centre. *)
let presence_bool =
  let m = Array.make_matrix tissue_X tissue_Y false in
  m.(tissue_X / 2).(tissue_Y / 2) <- true;
  m

(* AUTOMATICALLY COMPUTED presence matrix for cells. *)
let presence, cell_count =
  let m = Array.make_matrix (Array.length presence_bool) (Array.length presence_bool.(0)) (-1) in
  let c = ref 0 in
  for i = 0 to Array.length m - 1 do
    for j = 0 to Array.length m.(i) - 1 do
      if presence_bool.(i).(j) then
        begin
          m.(i).(j) <- diffusion_size + !c * variables_count;
          incr c
        end
      else ()
    done    
  done;
  m, !c

(* Initial state for cells + diffusion *)
let initial_state = 
  let stride = Array.length initial_state_per_cell in
  let tissue = Array.make (diffusion_size + cell_count * stride) initial_AHL_level in
  for i = 0 to tissue_Y - 1 do
    for j = 0 to tissue_X - 1 do
      let index = presence.(i).(j) in
      if index >= 0 then
        Array.blit initial_state_per_cell 0 tissue index variables_count
      else ()
    done
  done;
  tissue

let _ = Printf.printf "Model: %d cell(s), %d x %d space, %d variable(s).\n%!" cell_count tissue_X tissue_Y (Array.length initial_state)

(* An index to map rectangular coordinates to linear ones. *)
let index = indexing tissue_X tissue_Y

(* let _ = *)
(*   initial_state.(index.(10).(10)) <- 500.0 *)


(* Gradient temporary buffer. *)
let grad = Sundials.RealArray.make (tissue_X * tissue_Y) 0.0

let stop_time = float_of_string Sys.argv.(1)

(* Integration setup *)
let x0 = Nvector_serial.wrap (Sundials.RealArray.of_array initial_state)

let xt = Nvector_serial.wrap (Sundials.RealArray.of_array initial_state)

let tolerances = Cvode.SStolerances(1e-3, 1e-6)

let s = 
  let step = step tissue_X tissue_Y 0.1 index grad presence initial_parameters in
  Cvode.init Cvode.BDF (Cvode.Newton (Cvode.Dls.dense ())) tolerances step 0. x0

let _ = Cvode.set_stop_time s stop_time

(* Misc. *)

let compute_proteic_mass x k =
  compute_mM' x.{2} x.{3} x.{4} x.{5} x.{10} x.{11} x.{12} x.{13} k

let print_state = 
  let pattern =
    [| diffusion_size + var_s_i;
       diffusion_size + var_a;
       diffusion_size + var_r;
       diffusion_size + var_e_t;
       diffusion_size + var_e_m;
       diffusion_size + var_q;
       diffusion_size + var_m_r;
       diffusion_size + var_m_t;
       diffusion_size + var_m_m;
       diffusion_size + var_m_q;
       diffusion_size + var_c_r;
       diffusion_size + var_c_t;
       diffusion_size + var_c_m;
       diffusion_size + var_c_q;
       index.(tissue_X / 2).(tissue_Y / 2);
       diffusion_size + var_m_ahl;
       diffusion_size + var_c_ahl;
       diffusion_size + var_luxl;
       diffusion_size + var_m_luxl;
       diffusion_size + var_c_luxl;
       diffusion_size + var_luxr;
       diffusion_size + var_m_luxr;
       diffusion_size + var_c_luxr;
       diffusion_size + var_luxr_ahl
    |] in
  fprint_vec pattern "state"

let print_mass x  = fprint "mass" x

(* Graphics init. *)

let square_size = 10

let init () =
  let x = tissue_X * square_size in
  let y = tissue_Y * square_size in
  Graphics.open_graph "";
  Graphics.resize_window x y;
  Graphics.clear_graph ();
  Graphics.auto_synchronize false;
  Graphics.set_color Graphics.red

let _ = init ()

(* Integration *)
let _ =
  let state = Nvector.unwrap xt in
  display_diffusion tissue_X tissue_Y square_size index state;
  display_info 0;
  Graphics.synchronize ();
  for i = 1 to (int_of_float stop_time) do
    let (_, result) =
      try (Cvode.solve_normal s (float i) xt)
      with
      | exn -> raise exn
    in
    let _ =
      match result with
      | Cvode.Success ->
        () (* printf "success\n%!" *)
      | Cvode.RootsFound ->
        printf "roots found\n%!"
      | Cvode.StopTimeReached ->
        printf "stop time reached\n%!"
    in
    begin
      let state = Nvector.unwrap xt in
      if i mod 10 = 0 then
        begin
          display_diffusion tissue_X tissue_Y square_size index state;
          display_info i;
          Graphics.synchronize ()
        end;
      print_state (Some (Nvector.unwrap xt))
    end
    (* in *)
    (* let result_state = Nvector_serial.unwrap xt in *)
    (* print result_state *)
  done;
  print_state None


