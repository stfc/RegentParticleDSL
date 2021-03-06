-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

ast = require("regent/ast")
format = require("std/format")

compute_privileges = {}

local function return_true(node) return true end
local function return_false(node) return false end

local is_specialized_table = {
  [ast.specialized.stat.Reduce] = return_false,
  [ast.specialized.expr] = return_true,
  [ast.specialized.stat] = return_true,
  [ast.specialized.top] = return_true,
  [ast.specialized.Block] = return_true,
}

local is_spec_reduce_table = {
  [ast.specialized.expr.FieldAccess] = return_false,
  [ast.specialized.stat.Assignment] = return_false,
  [ast.specialized.expr] = return_true,
  [ast.specialized.stat] = return_true,
  [ast.specialized.top] = return_true,
  [ast.specialized.Block] = return_true,
}

local is_specialized_node = ast.make_single_dispatch(is_specialized_table, {}, return_false)()
local is_specialized_reduc_node = ast.make_single_dispatch(is_spec_reduce_table, {}, return_false)()

local function traverse_fieldaccess_postorder_two_region( node, sym1, sym2)
  local name = nil
  local symbol = 0
  ast.traverse_node_postorder(
    function(node)
        if(node:is(ast.specialized.expr.FieldAccess)) then
          if(node.value.value == sym1) then
            symbol = 1 
          elseif(node.value.value == sym2) then
            symbol = 2
          end
          if( name == nil) then
            name = node.field_name
          else
            name = name .. "." .. node.field_name
          end
        end
    end,
    node, is_specialized_node)
    if symbol == 0 then
      print("Error finding which symbol the FieldAccess equates to")
    end
    return name, symbol
end

local function traverse_fieldaccess_postorder_three_region( node, sym1, sym2, sym3)
  local name = nil
  local symbol = 0
  ast.traverse_node_postorder(
    function(node)
        if(node:is(ast.specialized.expr.FieldAccess)) then
          if(node.value.value == sym1) then
            symbol = 1
          elseif(node.value.value == sym2) then
            symbol = 2
          elseif(node.value.value == sym3) then
            symbol = 3
          end
          if( name == nil) then
            name = node.field_name
          else
            name = name .. "." .. node.field_name
          end
        end
    end,
    node, is_specialized_node)
    if symbol == 0 then
      print("Error finding which symbol the FieldAccess equates to")
    end
    return name, symbol

end

local function traverse_fieldaccess_postorder_one_region( node, sym1)
  local name = nil
  local symbol = 0
  ast.traverse_node_postorder(
    function(node)
        if(node:is(ast.specialized.expr.FieldAccess)) then
          if(node.value.value == sym1) then
            symbol = 1
          end
          if( name == nil) then
            name = node.field_name
          else
            name = name .. "." .. node.field_name
          end
        end
    end,
    node, is_specialized_node)
    if symbol == 0 then
      print("Error finding which symbol the FieldAccess equates to")
    end
    return name, symbol

end

local function traverse_specialized_postorder_two_region(fn, node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2, reduc_sym1, reduc_sym2)
  ast.traverse_node_postorder(
    function(node)
        fn(node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2, reduc_sym1, reduc_sym2)
    end,
    node, is_specialized_node)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2, reduc_sym1, reduc_sym2)
    end,
    node, is_specialized_reduc_node)
end


local function traverse_specialized_postorder_one_region(fn, node, sym, read_sym, write_sym, reduc_sym)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym, read_sym, write_sym, reduc_sym)
    end,
    node, is_specialized_node)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym, read_sym, write_sym, reduc_sym)
    end,
    node, is_specialized_reduc_node)
end

local function traverse_specialized_postorder_three_region(fn, node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3,
                                                           reduc_sym1, reduc_sym2, reduc_sym3)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3, reduc_sym1, reduc_sym2, reduc_sym3)
    end,
   node, is_specialized_node)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3, reduc_sym1, reduc_sym2, reduc_sym3)
    end,
   node, is_specialized_reduc_node)
end

local function three_region_privilege_map(node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3,
                                          reduc_sym1, reduc_sym2, reduc_sym3)
  --First case - we have an assignment. If we have an assignment and the lhs is a FieldAccess, then
  --this is a written to part of our region
  if(node:is(ast.specialized.stat.Assignment)) then
    if(node.lhs[1]:is(ast.specialized.expr.FieldAccess)) then
        local name, symbol = traverse_fieldaccess_postorder_three_region(node.lhs[1], sym1, sym2, sym3)
      if(symbol == 1) then
        write_sym1:insert(name)
      elseif(symbol == 2) then
        write_sym2:insert(name)
      elseif(symbol == 3) then
        write_sym3:insert(name)
      end
    end
  --Second case - For all field accesses we assume they are read from (though in some cases they are only
  --written to, I don't think that adding a read requirement should cause any issues)
  elseif(node:is(ast.specialized.expr.FieldAccess)) then
      if(node.value.value == sym1) then
        read_sym1:insert(node.field_name)
      elseif(node.value.value == sym2) then
        read_sym2:insert(node.field_name)
      elseif(node.value.value == sym3) then
        read_sym3:insert(node.field_name)
      end
  elseif(node:is(ast.specialized.stat.Reduce)) then
      local name, symbol = traverse_fieldaccess_postorder_three_region(node.lhs[1],sym1, sym2, sym3)
      if(symbol == 1) then
        reduc_sym1:insert({name, node.op})
      elseif(symbol == 2) then
        reduc_sym2:insert({name, node.op})
      elseif(symbol == 3) then
        reduc_sym3:insert({name, node.op})
      end
      --Recurse on RHS node to find any read field accesses.
      traverse_specialized_postorder_three_region(three_region_privilege_map, node.rhs, sym1, sym2, read_sym1,read_sym2,read_sym3, write_sym1, write_sym2, write_sym3,
                                                  reduc_sym1, reduc_sym2, reduc_sym3)
  end
end

local function two_region_privilege_map(node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2, reduc_sym1, reduc_sym2)
  --First case - we have an assignment. If we have an assignment and the lhs is a FieldAccess, then
  --this is a written to part of our region
  if(node:is(ast.specialized.stat.Assignment)) then
    if(node.lhs[1]:is(ast.specialized.expr.FieldAccess)) then
        local name, symbol = traverse_fieldaccess_postorder_two_region(node.lhs[1],sym1, sym2)
      if(symbol == 1) then
        write_sym1:insert(name)
      elseif(symbol == 2) then
        write_sym2:insert(name)
      end
    else
        --For array accesses, we recurse until we find a non-array access and then check if its a FieldAccess
        if node.lhs[1]:is(ast.specialized.expr.IndexAccess) then
            local z = node.lhs[1]
            while z.value:is(ast.specialized.expr.IndexAccess) do
                z = z.value
            end
            if z.value:is(ast.specialized.expr.FieldAccess) then
                local name, symbol = traverse_fieldaccess_postorder_two_region(z.value,sym1, sym2)
                if(symbol == 1) then
                    write_sym1:insert(name)
                elseif(symbol == 2) then
                    write_sym2:insert(name)
                end
            end
        end
    end
  --Second case - For all field accesses we assume they are read from (though in some cases they are only
  --written to, I don't think that adding a read requirement should cause any issues)
  elseif(node:is(ast.specialized.expr.FieldAccess)) then
      if(node.value.value == sym1) then
        read_sym1:insert(node.field_name)
      elseif(node.value.value == sym2) then
        read_sym2:insert(node.field_name)
      end
  elseif(node:is(ast.specialized.stat.Reduce)) then
      local name, symbol = traverse_fieldaccess_postorder_two_region(node.lhs[1],sym1, sym2)
      if(symbol == 1) then
        reduc_sym1:insert({name, node.op})
      elseif(symbol == 2) then
        reduc_sym2:insert({name, node.op})
      end
      --Recurse on RHS node to find any read field accesses.
      traverse_specialized_postorder_two_region(two_region_privilege_map, node.rhs, sym1, sym2, read_sym1,read_sym2, write_sym1, write_sym2,
                                                reduc_sym1, reduc_sym2)
  end
end

local function one_region_privilege_map(node, sym, read_sym, write_sym, reduc_sym)
  --First case - we have an assignment. If we have an assignment and the lhs is a FieldAccess, then
  --this is a written to part of our region
  if(node:is(ast.specialized.stat.Assignment)) then
    if(node.lhs[1]:is(ast.specialized.expr.FieldAccess)) then
      local name, symbol = traverse_fieldaccess_postorder_one_region(node.lhs[1],sym)
      if(symbol == 1) then
        write_sym:insert(name)
      end
    end
  --Second case - For all field accesses we assume they are read from (though in some cases they are only
  --written to, I don't think that adding a read requirement should cause any issues)
  elseif(node:is(ast.specialized.expr.FieldAccess)) then
      if(node.value.value == sym) then
        read_sym:insert(node.field_name)
      end
  elseif(node:is(ast.specialized.stat.Reduce)) then
      local name, symbol = traverse_fieldaccess_postorder_one_region(node.lhs[1],sym)
      if(symbol == 1) then
        reduc_sym:insert({name, node.op})
      end
      --Recurse on RHS node to find any read field accesses.
      traverse_specialized_postorder_one_region(one_region_privilege_map, node.rhs, sym, read_sym, write_sym, reduc_sym)
  end
end



function compute_privileges.three_region_privileges(kernel_name)
  local part1 = regentlib.newsymbol("part1")
  local part2 = regentlib.newsymbol("part2")
  local r2 = regentlib.newsymbol("r2")
  local config = regentlib.newsymbol("config")
  local read1 = terralib.newlist()
  local read2 = terralib.newlist()
  local read3 = terralib.newlist()
  local write1 = terralib.newlist()
  local write2 = terralib.newlist()
  local write3 = terralib.newlist()
  local reduc1 = terralib.newlist()
  local reduc2 = terralib.newlist()
  local reduc3 = terralib.newlist()
  local kernel = kernel_name(part1, part2, r2, config)
  traverse_specialized_postorder_three_region( three_region_privilege_map, kernel.ast, part1, part2, config, read1, read2, read3, write1, write2, write3,
                                                reduc1, reduc2, reduc3 )
  local hash = {}
  local temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(read1) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  read1 = temp

  hash = {}
  temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(read2) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  read2 = temp

  hash = {}
  temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(read3) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  read3 = temp

  hash = {}
  temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(write1) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  write1 = temp

  hash = {}
  temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(write2) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  write2 = temp

  hash = {}
  temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(write3) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  write3 = temp

  temp = terralib.newlist()
  for _, v in pairs(reduc1) do
    local exists = false
    for _, v2 in pairs(temp) do
      if v[1] == v2[1] and v[2] == v2[2] then
        exists = true
      end
    end
    if( not exists) then
      temp:insert(v)
    end
  end
  reduc1 = temp
  temp = terralib.newlist()
  for _, v in pairs(reduc2) do
    local exists = false
    for _, v2 in pairs(temp) do
      if v[1] == v2[1] and v[2] == v2[2] then
        exists = true
      end
    end
    if( not exists) then
      temp:insert(v)
    end
  end
  reduc2 = temp
  temp = terralib.newlist()
  for _, v in pairs(reduc3) do
    local exists = false
    for _, v2 in pairs(temp) do
      if v[1] == v2[1] and v[2] == v2[2] then
        exists = true
      end
    end
    if( not exists) then
      temp:insert(v)
    end
  end
  reduc3 = temp

  --Check for collisions between reduction and read/writes
  hash = {}
  for _,v in pairs(reduc1) do
    hash[v[1]] = true
  end
  for _,v in pairs(read1) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a read and reduction state. This may impact performance.")
    end
  end
  for _,v in pairs(write1) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a write and reduction state. This may impact performance.")
    end
  end

  --Check for collisions between reduction and read/writes
  hash = {}
  for _,v in pairs(reduc2) do
    hash[v[1]] = true
  end
  for _,v in pairs(read2) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a read and reduction state. This may impact performance.")
    end
  end
  for _,v in pairs(write2) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a write and reduction state. This may impact performance.")
    end
  end

  --Check for collisions between reduction and read/writes
  hash = {}
  for _,v in pairs(reduc3) do
    hash[v[1]] = true
  end
  for _,v in pairs(read3) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a read and reduction state. This may impact performance.")
    end
  end
  for _,v in pairs(write3) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a write and reduction state. This may impact performance.")
    end
  end


  return read1, read2, read3, write1, write2, write3, reduc1, reduc2, reduc3
end

function compute_privileges.one_region_privileges(kernel_name)
  local part1 = regentlib.newsymbol("part")
  local z = kernel_name(part1)
  local read = terralib.newlist()
  local write = terralib.newlist()
  local reduc = terralib.newlist()
  traverse_specialized_postorder_one_region( one_region_privilege_map, z.ast, part1, read, write, reduc)
  local hash = {}
  local temp = terralib.newlist()
  --Remove repeats before returning
  for _,v in pairs(read) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  read = temp
  hash = {}
  temp = terralib.newlist()
  for _,v in pairs(write) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  write = temp
  temp = terralib.newlist()
  for _, v in pairs(reduc) do
    local exists = false
    for _, v2 in pairs(temp) do
      if v[1] == v2[1] and v[2] == v2[2] then
        exists = true
      end
    end
    if( not exists) then
      temp:insert(v)
    end
  end
  reduc = temp

  --Check for collisions between reduction and read/writes
  hash = {}
  for _,v in pairs(reduc) do
    hash[v[1]] = true
  end
  for _,v in pairs(read) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a read and reduction state. This may impact performance.")
    end
  end
  for _,v in pairs(write) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a write and reduction state. This may impact performance.")
    end
  end

  return read,write, reduc
end

function compute_privileges.two_region_privileges(kernel_name)
  local part1 = regentlib.newsymbol("part1")
  local part2 = regentlib.newsymbol("part2")
  local r2 = regentlib.newsymbol("r2")
--  print(kernel_name(part1, part2, r2))
--  print(kernel_name)
  local z = kernel_name(part1, part2, r2)
--  print(z.ast)
  local read_p1 = terralib.newlist()
  local read_p2 = terralib.newlist()
  local write_p1 = terralib.newlist()
  local write_p2 = terralib.newlist()
  local reduc_p1 = terralib.newlist()
  local reduc_p2 = terralib.newlist()
  traverse_specialized_postorder_two_region(  two_region_privilege_map, z.ast, part1, part2, read_p1, read_p2, write_p1, write_p2, reduc_p1, reduc_p2)
  --Remove repeats before returning
  local hash = {}
  local temp = terralib.newlist()
  for _,v in pairs(read_p1) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  read_p1 = temp
  hash = {}
  temp = terralib.newlist()
  for _,v in pairs(read_p2) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  read_p2 = temp
  hash = {}
  temp = terralib.newlist()
  for _,v in pairs(write_p1) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  write_p1 = temp
  hash = {}
  temp = terralib.newlist()
  for _,v in pairs(write_p2) do
    if( not hash[v]) then
      temp:insert(v)
      hash[v] = true
    end
  end
  write_p2 = temp

  temp = terralib.newlist()
  for _, v in pairs(reduc_p1) do
    local exists = false
    for _, v2 in pairs(temp) do
      if v[1] == v2[1] and v[2] == v2[2] then
        exists = true
      end
    end
    if( not exists) then
      temp:insert(v)
    end
  end
  reduc_p1 = temp
  temp = terralib.newlist()
  for _, v in pairs(reduc_p2) do
    local exists = false
    for _, v2 in pairs(temp) do
      if v[1] == v2[1] and v[2] == v2[2] then
        exists = true
      end
    end
    if( not exists) then
      temp:insert(v)
    end
  end
  reduc_p2 = temp

  --Check for collisions between reduction and read/writes
  hash = {}
  for _,v in pairs(reduc_p1) do
    hash[v[1]] = true
  end
  for _,v in pairs(read_p1) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a read and reduction state. This may impact performance.")
    end
  end
  for _,v in pairs(write_p1) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a write and reduction state. This may impact performance.")
    end
  end

  hash = {}
  for _,v in pairs(reduc_p2) do
    hash[v[1]] = true
  end
  for _,v in pairs(read_p2) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a read and reduction state. This may impact performance.")
    end
  end
  for _,v in pairs(write_p2) do
    if(hash[v]) then
      print("WARNING: Access to field", v, "appears in both a write and reduction state. This may impact performance.")
    end
  end
  return read_p1, read_p2, write_p1, write_p2, reduc_p1, reduc_p2
end




return compute_privileges
