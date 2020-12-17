import "regent"

get_args = require("src/utils/read_args")

--Choose options for WCSPH interactions.
local interactions = {}
local interactions_string = {}
local interactions_LAMINAR_SPS = 1
local interactions_DEFAULT = interactions_LAMINAR_SPS
interactions[interactions_DEFAULT] = "src/interactions/WC_SPH/force"
interactions_string["LAM_SPS"] = interactions_LAMINAR_SPS
interactions_string["DEFAULT"] = interactions_DEFAULT



local shifting = {}
local shifting_string = {}
local shifting_NONE = 1
local shifting_DEFAULT = shifting_NONE
shifting[shifting_NONE] = "src/interactions/WC_SPH/no_shifting"
shifting_string["NONE"] = shifting_NONE
shifting_string["DEFAULT"] = shifting_NONE


local diffusion = {}
local diffusion_string = {}
local diffusion_NONE = 1
local diffusion_DEFAULT = diffusion_NONE
diffusion[diffusion_DEFAULT] = "src/interactions/WC_SPH/no_density_diffusion"
diffusion_string["NONE"] = diffusion_NONE
diffusion_string["DEFAULT"] = diffusion_DEFAULT

local boundary = {}
local boundary_string = {}
local boundary_SIMPLE = 1
local boundary_DEFAULT = boundary_SIMPLE
boundary[boundary_SIMPLE] = "src/interactions/WC_SPH/simple_boundary"
boundary_string["SIMPLE"] = boundary_SIMPLE
boundary_string["DEFAULT"] = boundary_DEFAULT

local interactions_chosen = interactions[interactions_DEFAULT]
local shifting_chosen = shifting[shifting_DEFAULT]
local diffusion_chosen = diffusion[diffusion_DEFAULT]
local boundary_chosen = boundary[boundary_DEFAULT]

--Read in some parameters from the command line.
local function select_SPH_method()
  local arg_ints = get_args.get_optional_arg("-wcsph:int")
  local arg_shift = get_args.get_optional_arg("-wcsph:shifting")
  local arg_diff = get_args.get_optional_arg("-wcsph:diffusion")
  local arg_bound = get_args.get_optional_arg("-wcsph:boundary")
  --Check the interaction type
  if arg_ints ~= nil then
    --Check if its a known interaction type
    local i = interactions_string[arg_ints]
    if i ~= nil then
      local x = interactions[i]
      if x ~= nil then
        interactions_chosen = x
      end
    end
  end

  --Check the shifting type
  if arg_shift ~= nil then
    local i = shifting_string[arg_shift]
    if i ~= nil then
      local x = shifting[i]
      if x ~= nil then
        shifting_chosen = x
      end
    end
  end

  --Check the diffusion type
  if arg_diff ~= nil then
    local i = diffusion_string[arg_diff]
    if i ~= nil then
      local x = diffusion[i]
      if x ~= nil then
        diffusion_chosen = x
      end
    end
  end

  --Check the boundary type
  if arg_bound ~= nil then
    local i = boundary_string[arg_bound]
    if i ~= nil then
      local x = boundary[i]
      if x ~= nil then
        boundary_chosen = x
      end
    end
  end
end


select_SPH_method()

require(shifting_chosen)
require(diffusion_chosen)
require(boundary_chosen)
require(interactions_chosen)
