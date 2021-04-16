import "regent"

--dl_meso's config type contains a variety of extra things needed for dl_meso simulations

CONST_max_species = 50
CONST_max_mxprm = 7
CONST_stksize = 13
CONST_statsize = 52
CONST_maxstk = 101

CONST_randsize = 624

--FIXME: Remove this to be part of the main config import in RegentParticleDSL
require("src/config/space")

fspace config_type{
  space : space_config_type,
--  neighbour_config : neighbour_config_type,
--  timing_config : timing_config_type,

--Used for scan/read field
  cutoff : double,
--  srfzcut : double,
  mxprm : int,
  ltabpot : bool,
  lrcut : bool,
  nspe : int,
  npot : int,
  namtmp : int8[CONST_max_species][9],
  masstmp : double[CONST_max_species],
  chgetmp : double[CONST_max_species],
  eunit : double,
  gamma : double[CONST_max_mxprm],
  vvv : double[CONST_max_mxprm][CONST_max_species * (CONST_max_species+1) / 2],
  ktype : int[CONST_max_species],

  --Scan CONTROL fields
  l_exist : bool,
  l_safe : bool,
  l_scr : bool,
  l_temp : bool,
  l_time : bool,
  l_conf : bool,
  l_init : bool,
  l_rest : bool,
  temp : double,
  tstep : double,
  rtstep : double,
  tstepsq : double,

  --Read CONTROL fields
  tclose : double,
  itype : int,
  btype : int,
  nseql : int,
  ldyn : bool,
  nsbpo : int,
  nsbts : int, 
  ltemp : bool,
  nstk : int,
  iscorr : int,
  lcorr : bool,
  nrun : int, 
  straj : int,
  ntraj : int,
  keytrj : int,
  ltraj : bool,
  volm : double,

  outsel : double,

  --Step counter
  nstep : int,

  --OUTPUT fields - for now located here.
  stpte : double,
  stppe : double,
  stpee : double,
  stpbe : double,
  stpae : double,
  stpde : double,
  stpse : double,
  stpvir : double,
  stptke : double,
  stpprs : double,
  stpvlm : double,
  stpzts : double,
  stpttp: double,
  stptpx : double,
  stptpy : double,
  stptpz : double,
  stpbdl : double,
  stpang : double,
  stpdhd : double,
  rav : double[CONST_stksize],
  ave : double[CONST_statsize],
  flc : double[CONST_statsize],
  zum : double[12],

  --CORREL field.
  CORREL_newjob : bool,
  stress : double[36],

  --HISTORY fields
  kres : int,
  nusyst: int, --number of unbonded beads
  nsyst : int, --number of total beads


  --Header size in HISTORY file
  headersize : int64,
  numframe : int32,
  filesize : int64,
  markerpos : int64,

  --Displacement of Lees-Edwards shearing boundary
  shrdx : double,
  shrdy : double,
  shrdz : double,


  --Simulation name
  text : int8[80],

  --Barostat properties and statistical accumulators
    --Barostat piston velocity
    upx : double,
    upy : double,
    upz : double,
    --Barostat piston force
    fpx : double,
    fpy : double,
    fpz : double,

    --Current number of timesteps for statistical sampling
    nav :  int32,
  
    --Statistics arrays
    stkpe : double[CONST_maxstk],
    stkee : double[CONST_maxstk],
    stkse : double[CONST_maxstk],
    stkde : double[CONST_maxstk],
    stkae : double[CONST_maxstk],
    stkbe : double[CONST_maxstk],
    stkvir : double[CONST_maxstk],
    stkvlm : double[CONST_maxstk],
    stkzts : double[CONST_maxstk],
    stktkex : double[CONST_maxstk],
    stktkey : double[CONST_maxstk],
    stktkez : double[CONST_maxstk],

    --Random state
    mt : int32[CONST_randsize],
}
