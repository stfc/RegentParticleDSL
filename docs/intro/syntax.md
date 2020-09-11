# Regent Language

This section of the documentation is dedicated to features of the Regent programming language that differ
from more commonly used languages, such as Python or C.

## Headers

All files used with RegentParticleDSL must begin with the same first line:
```
    import "regent"
```

This tells the regent compiler to recognize the file correctly.

To include code or declarations from other files, the language uses the require code. For example:
```
    require("path/to/file")
```

will import the file at `path/to/file.rg` from the current directory. I think absolute paths can
also be used, however relative paths using dots (e.g. `../path`) are not possible.


## Comments

Comments in Regent are inline comments, and start with `--`. 

## Local and global variables

Variables and functions declared in Regent are usually global in scope when outside of functions.
This can be avoided by using the `local` flag, e.g:
```
    local variable_name = 1.0
    
    local function func()
    --Code goes here
    end
```

## Regent or Lua

Most code users write in RegentParticleDSL is Lua, with only code inside particle structures, kernels and
the main program having Regent's specific syntax. Regent variables are declared as:
```
    var int_variable : int = 1
    var float_variable : float = 1.0
```

The other features specific to Regent are discussed later parts of the documentation.

## Terra
As well as Regent and Lua code, RegentParticleDSL also supports Terra code inside files: [terralang.org](http://terralang.org/).
