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
  [ast.specialized.expr] = return_true,
  [ast.specialized.stat] = return_true,
  [ast.specialized.top] = return_true,
  [ast.specialized.Block] = return_true,
}

local is_specialized_node = ast.make_single_dispatch(is_specialized_table, {}, return_false)()

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


local function three_region_privilege_map(node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3)
  --First case - we have an assignment. If we have an assignment and the lhs is a FieldAccess, then
  --this is a written to part of our region
  if(node:is(ast.specialized.stat.Assignment)) then
    if(node.lhs[1]:is(ast.specialized.expr.FieldAccess)) then
        local name, symbol = traverse_fieldaccess_postorder_three_region(node.lhs[1], sym1, sym2)
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
  end
end

local function two_region_privilege_map(node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2)
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
    end
  --Second case - For all field accesses we assume they are read from (though in some cases they are only
  --written to, I don't think that adding a read requirement should cause any issues)
  elseif(node:is(ast.specialized.expr.FieldAccess)) then
      if(node.value.value == sym1) then
        read_sym1:insert(node.field_name)
      elseif(node.value.value == sym2) then
        read_sym2:insert(node.field_name)
      end
  end
end

local function one_region_privilege_map(node, sym, read_sym, write_sym)
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
  end
end


local function traverse_specialized_postorder_two_region(fn, node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2)
  ast.traverse_node_postorder(
    function(node)
        fn(node, sym1, sym2, read_sym1, read_sym2, write_sym1, write_sym2)
    end,
    node, is_specialized_node)
end


local function traverse_specialized_postorder_one_region(fn, node, sym, read_sym, write_sym)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym, read_sym, write_sym)
    end,
    node, is_specialized_node)
end

local function traverse_specialized_postorder_three_region(fn, node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3)
  ast.traverse_node_postorder(
    function(node)
      fn(node, sym1, sym2, sym3, read_sym1, read_sym2, read_sym3, write_sym1, write_sym2, write_sym3)
    end,
   node, is_specialized_node)
end

function compute_privileges.three_region_privileges(kernel_name)
  local r1 = regentlib.newsymbol("region1")
  local r2 = regentlib.newsymbol("region2")
  local r3 = regentlib.newsymbol("region3")
  local read1 = terralib.newlist()
  local read2 = terralib.newlist()
  local read3 = terralib.newlist()
  local write1 = terralib.newlist()
  local write2 = terralib.newlist()
  local write3 = terralib.newlist()
  local kernel = kernel_name(r1, r2, r3)
  traverse_specialized_postorder_three_region( three_region_privilege_map, kernel.ast, r1, r2, r3, read1, read2, read3, write1, write2, write3 )
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


  return read1, read2, read3, write1, write2, write3
end

function compute_privileges.one_region_privileges(kernel_name)
  local part1 = regentlib.newsymbol("part")
  local z = kernel_name(part1)
  local read = terralib.newlist()
  local write = terralib.newlist()
  traverse_specialized_postorder_one_region( one_region_privilege_map, z.ast, part1, read, write)
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
  return read,write
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
  traverse_specialized_postorder_two_region(  two_region_privilege_map, z.ast, part1, part2, read_p1, read_p2, write_p1, write_p2)
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
  return read_p1, read_p2, write_p1, write_p2
end




return compute_privileges
