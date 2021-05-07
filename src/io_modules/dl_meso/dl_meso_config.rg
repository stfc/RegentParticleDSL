import "regent"

--dl_meso's config type contains a variety of extra things needed for dl_meso simulations

CONST_max_species = 50
CONST_max_potential = CONST_max_species * (CONST_max_species + 1) /2
CONST_max_mxprm = 7
CONST_stksize = 13
CONST_statsize = 52
CONST_maxstk = 101

CONST_randsize = 624

--FIXME: Remove this to be part of the main config import in RegentParticleDSL
require("src/config/space")

fspace force_mdvv_type{
    strscxx : double,
    strscxy : double,
    strscxz : double,
    strscyy : double,
    strscyz : double,
    strsczz : double,
    strsdxx : double,
    strsdxy : double,
    strsdxz : double,
    strsdyy : double,
    strsdyz : double,
    strsdzz : double,
    strsrxx : double,
    strsrxy : double,
    strsrxz : double,
    strsryy : double,
    strsryz : double,
    strsrzz : double,
}

fspace config_type{
  space : space_config_type,
  neighbour_config : neighbour_config_type,
--  timing_config : timing_config_type,

--Used for scan/read field
  cutoff : double,
  rtcut : double,
  rtct2 : double,
  rrtcut : double,
  rct2 : double,
  rrct2 : double,
  gamma : double[CONST_max_potential],
  sigma : double[CONST_max_potential],
--  srfzcut : double,
  mxprm : int32,
  ltabpot : bool,
  lrcut : bool,
  lvarfc : bool,
  nspe : int32,
  npot : int32,
  namspe : int8[CONST_max_species][9],
  masstmp : double[CONST_max_species],
  chgetmp : double[CONST_max_species],
  nspec : int[CONST_max_species],
  nusystcell : int32,
  nsystcell : int32,
  eunit : double,
--  gamma : double[CONST_max_mxprm],
  vvv : double[CONST_max_mxprm][CONST_max_potential],
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

--Time fields
  timfrc : double,
  timstp : double,

  l_config : bool,

  --Timing fields
  time : int64,
  tzero : int64,

  --Read CONTROL fields
  tclose : double,
  itype : int32,
  btype : int32,
  ldyn : bool,
  nsbpo : int32,
  nsbts : int32, 
  ltemp : bool,
  nstk : int32,
  iscorr : int32,
  lcorr : bool,
  nrun : int32, 
  straj : int32,
  ntraj : int32,
  keytrj : int32,
  ltraj : bool,
  volm : double,

  outsel : double,

  --Step counter
  nstep : int32,
  nseql : int32, --Number of timesteps for system equilibration
  ndump : int32, --Frequency for creating simulation restart data

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
  stptkex : double,
  stptkey : double,
  stptkez : double,
  stpbdl : double,
  stpbdmx : double,
  stpbdmn : double,
  stpang : double,
  stpdhd : double,
  rav : double[CONST_stksize],
  ave : double[CONST_statsize],
  flc : double[CONST_statsize],
  zum : double[12],

  --CORREL field.
  CORREL_newjob : bool,

  --HISTORY fields
  kres : int32,
  nusyst: int32, --number of unbonded beads
  nsyst : int32, --number of total beads
  nfsyst : int32, --?


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
  text : int8[81],

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

    --barostat parameters
    abaro : double,
    psmass : double,
    rpsmass : double,
    sigmalang : double,
    fkt : double,
    bbaro : double,

    --Statisical properties used during interactions
    pe : double, --potential energy
    vir : double, --virial
    ivrl : double[3], --Instantaenous virial
    be : double, --bond energy
    ae : double, --angle energy
    de : double, --dihedral energy
    ee : double, --electrostatic energy
    se : double, --surface energy
    bdlng : double, --summed bond length
    bdlmin : double, --Minimum bond length
    bdlmax : double, --Maximum bond length
    bdang : double, --Summed bond angles
    bddhd : double, -- Summed bond dihedrals
    tke : double[3], --Kinetic energy separated into x-, y- and z-components
    stress : double[36], -- Stress tensor separated into conservative, dissipative, random and
                        -- kinetic contributions

    mdvv_type : force_mdvv_type,
}
