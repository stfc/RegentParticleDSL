------------------------------------------------------------- 
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

local io_utils = {}
local c_string = terralib.includec("string.h")
local c_ctype = terralib.includec("ctype.h")
local c_stdlib = terralib.includec("stdlib.h")
local c_stdio = terralib.includec("stdio.h")

--Max word length used for parsing - one higher than fortran as we manually need space for '\0'
local mxword = 40
--Gets the nth word from a line of delimited text.
--@param text : Line of text to be read (int8acter array)
--@param n : The number of the word to be returned
--@return &int8: The word to be returned. NB: This needs to be manually free'd to avoid memory leakage!
terra get_word( text : &int8, n : int) : &int8
  var cur_word = 0
  var a = 0
  var b = 0 
  var position = 0
  var length = c_string.strlen(text)

  var word : &int8 = [&int8](regentlib.c.malloc( mxword ))
  while (cur_word < n and position < length) do
    var u = text[position]
    --Check for space, comma or tab
    if (u == int8(32) or u == int8(44) or u == int8(9)) then
      	if b > a then
        	cur_word = cur_word + 1
      	end
      	if cur_word < n then
        	a = position + 1
      	end
    else
 		b = position + 1
    end
    position = position  + 1 
  end

  if cur_word == n then
        var ii = 0
        for i = a,b do
            word[ii] = text[i]
            ii = ii + 1
        end
        word[ii] = int8(0)
  elseif position == length and cur_word == n-1 then
        var ii = 0
        for i = a, length+1 do
            word[ii] = text[i]
            ii = ii + 1
        end
        word[ii] = int8(0)
  end

  var ii = 0
  while word[ii] ~= int8(0) do
    if word[ii] == int8(10) then
        word[ii] = int8(0)
    end
    ii = ii+1
  end
  return word
end

io_utils.get_word = get_word

terra io_utils.to_lowercase(text : &&int8)

  var i = 0
  while (@text)[i] ~= int8(0) do
    (@text)[i] = c_ctype.tolower((@text)[i])
    i = i + 1
  end

end

terra io_utils.get_double(text : &int8, n : int ) : double

  var word : &int8 = get_word(text, n)
  var value : double = c_stdlib.strtod(word, [&&int8](0))
  regentlib.c.free(word)
  return value
end

terra io_utils.ANY( logical : &&bool, s1: int, s2: int, value : bool) : bool

  var i : int
  var j : int
  var result : bool = true
  for i=0, s1 do
    for j=0, s2 do
        result = result and (logical[i][j] == value)
    end
  end
  return result
end


terra io_utils.get_file_size( filename : rawstring ) : int64
  var f = c_stdio.fopen(filename, 'rb')
  c_stdio.fseek(f, 0,  c_stdio.SEEK_END)
  var size : int32 = c_stdio.ftell(f)
  c_stdio.fclose(f)
  return size
end

return io_utils
