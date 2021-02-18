import "regent"

local format = require("std/format")

if DSL_settings.TIMING then
  fspace timing_config_type{
    user_task_time : int64,
    neighbour_search_time: int64,
    other_time: int64
  }
else
  fspace timing_config_type{

  }
end

function DSL_report_timings()
  local timings_quote = rquote

  end
  if DSL_settings.TIMING then
    timings_quote = rquote
      format.println("Time in user tasks: {}s", double([variables.config][0].timing_config.user_task_time) / 1000000.0)
      format.println("Time in neighbour search overhead: {}s", double([variables.config][0].timing_config.neighbour_search_time) / 1000000.0)
      format.println("Time in other areas: {}s", double([variables.config][0].timing_config.other_time) / 1000000.0)
    end
  end
  return timings_quote
end
