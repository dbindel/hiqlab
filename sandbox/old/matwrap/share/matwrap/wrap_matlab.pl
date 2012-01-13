#
# Language module for the wrapper generator for matlab.
#
# MATLAB is a little difficult to handle because it can only define a single
# C function per mex file.  This means that the C function must take an
# additional argument of which sub-function to call.  We provide MATLAB
# functions (each in their own addition file--argh!!!) which supply the
# extra argument.
#
# Copyright (c) 1997 Gary R. Holt.  This is distributed under the terms of the 
# perl artistic license (http://language.perl.com/misc/Artistic.html).
#
# This code has slightly funny stylistic conventions, largely so I can get
# emacs perl mode to get the indentation correct.  In particular, the use
# of here documents and multi-line double quoted strings inside parentheses is
# probably a little weird.
#

package matlab;

$outdir = '.';			# Default directory for the other files that
				# we have to create.
@function_names = ();		# This table translates function names and
				# their assigned IDs.  The extra numerical
				# argument to the mexFunction() indicates
				# which function the user actually wanted to
				# call.

$max_dimensions = 4;		# Maximum number of dimensions.  (Although
				# matlab doesn't have a fixed number of
				# dimensions, we only can handle up to
				# 4 dimensional arrays.  If you need more,
				# just up this number; no other modifications
				# are necessary.)

*pointer_type_code = *main::pointer_type_code;
				# Copy the definition of pointer_type_code
				# into this module.

$use_mxCalloc = 0;		# By default, use alloca for allocating memory.
#
# arg_pass(\%function_def, $argname)
#
# A C or C++ expression used to pass the argument to another function
# which does not know anything about the type of the argument.  For
# example, in the MATLAB module this function returns an expression for
# the mxArray type for a given argument.
#
sub arg_pass {
  my ($faa, $argname) = @_;	# Name the arguments.
  $faa->{args}{$argname}{mat_expr}; # Return the octave expression.
}

#
# arg_declare("arg_name_in_arglist")
# 
# This returns a C/C++ declaration appropriate for the argument passed
# using arg_pass.  For example, in the MATLAB module this function returns
# "mxArray *arg_name_in_arglist".
#
sub arg_declare {
  "const mxArray *$_[0]";
}

#
# declare_const("constant name", "class name", "<type>", "doc str")
#
# Output routines to make a given constant value accessible from the 
# interpreter.
# If "class name" is blank, this is a global constant.
#
sub declare_const {
				# This is currently not supported.
}

#
# error_dimension(\%function_def, $argname)
#
# A C statement (including the final semicolon, if not surrounded by braces)
# which indicates that an error has occured because the dimension of argument
# $argname was wrong.
#
sub error_dimension {
  my ($faa, $argname) = @_;	# Name the arguments.
  "mexErrMsgTxt(\"Error in dimension of argument $argname\");";
}

#
# finish()
#
# Called after all functions have been wrapped, to close the output file and do
# whatever other cleanup is necessary.
#
# For MATLAB, we now have to dispatch to the appropriate wrapper function
# that we have created.
#
sub finish {

  print OUTFILE "

/*
 * The main dispatch function.  This function calls the appropriate wrapper
 * based on the value of the first argument.
 */

#ifdef __cplusplus
extern \"C\" {
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs == 0 || !mxIsDouble(prhs[0]) ||
      mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1)
    mexErrMsgTxt(\"Illegal call to mex functions in $outfile\");
  first_arg = prhs[0];
  switch ((int)(mxGetScalar(prhs[0]))) {
";
				# Output the beginning of the definition.
  foreach (0 .. @function_names-1) {
    print "  case $_: $function_names[$_](nlhs, plhs, nrhs, prhs); break;\n";
  }
				# Dispatch appropriately.
  print "
  default: mexErrMsgTxt(\"Illegal function ID parameter\");
  }
}

#ifdef __cplusplus
}
#endif
";
  
  close(OUTFILE);		# Done with this file.
}

#
# Begin the definition of a function.  Arguments:
# 1) The %function_def array.
#
sub function_start {
  my ($faa) = @_;		# Access the argument.

  my $fname = $faa->{script_name} ||
    ($faa->{class} ? $faa->{class} . "_" : "") . $faa->{name};
				# Compute the name from MATLAB.
  $faa->{matlab_name} = $fname;	# Remember it.
#
# Associate each argument to the function with an index:
#
  my $argno = 1;		# Start numbering input/modify arguments at
				# 1, since the first argument is the function
				# ID.
  foreach (@{$faa->{inputs}}, @{$faa->{modifies}}) {
    my $arg = $faa->{args}{$_};	# Point to the argument array.
    $arg->{mat_expr} = "prhs[" . $argno++ . "]";
    if ($arg->{basic_type} eq 'char *' || $arg->{basic_type} eq 'm_Object') {
      if (@{$arg->{dimension}} != 0) { # A vector of strings?
	die("wrap_matlab: vectors of strings are not currently supported\n");
      }
      $arg->{vectorize} = 0;	# We don't support vectorization of char *.
      $faa->{vectorize} = 0 if $arg->{source} eq 'modify';
				# If it's a modify variable, we can't vectorize
				# the function at all.
    }
  }

  $argno = 0;			# Output arguments are numbered starting at 0.
  foreach (@{$faa->{outputs}}) {
    my $arg = $faa->{args}{$_};	# Point to the argument array.
    $arg->{mat_expr} = "plhs[" . $argno++ . "]";
    if ($arg->{basic_type} eq 'char *' || $arg->{basic_type} eq 'm_Object') {
      if (@{$arg->{dimension}} != 0) { # A vector of characters?
	die("wrap_matlab: vectors of strings are not currently supported\n");
      }
      $faa->{vectorize} = 0;	# This function cannot be vectorized if it has
				# string output.
    }
  }


  my $plhs = (@{$faa->{outputs}} ? 'plhs' : '');
  my $prhs = (@{$faa->{inputs}} + @{$faa->{modifies}} ? 'prhs' : '');
  ("void _wrap_$fname(int nlhs, mxArray **$plhs, int nrhs, const mxArray **$prhs)\n" .
   "{\n" .
   (@{$faa->{outputs}} == 1 ?
    "  if (nlhs > 1 ||\n" :
    "  if (nlhs != " . @{$faa->{outputs}} . " ||\n") .
   "      nrhs != " . (@{$faa->{inputs}} + @{$faa->{modifies}}) . "+1)\n" .
				# +1 is because the first argument to the
				# function is actually the function code.
   "    mexErrMsgTxt(\"Wrong number of arguments to function $fname\");\n");
				# Start the function definition.
}

#
# function_end(\%function_def)
#
# Return a string which finishes off the definition of a wrapper.
#
sub function_end {
  my ($faa) = @_;		# Access the argument.

  my $fname = $faa->{matlab_name}; # Get matlab's name for the function.
  push(@function_names, "_wrap_$fname"); # Register it in our list of functions
  my $doc = $faa->{comment};

#
# We also have to generate an auxiliary .m file which contains the help text
# and a wrapper for the function:
#
  local (*MFILE);		# Make a local file handle.
  open(MFILE, "> $outdir/$fname.m") || # Create the file.
    die("can't open file $outdir/$fname.m--$!\n");

  my $outargs = join(", ", @{$faa->{outputs}});	# Make list of output arguments.
  @{$faa->{outputs}} > 1 and $outargs = "[$outargs]";
				# If there's more than one, put it in
				# brackets.
  $outargs ne '' and $outargs .= " = ";	# Add an equals sign too.
  my @inargs = @{$faa->{inputs}}; # Make a copy of the input arguments.
  foreach (@inargs) { s/^_+//; } # Strip off leading underscores.  For some
				# strange reason matlab doesn't like these.

  my $inargs = join(", ", @inargs, @{$faa->{modifies}});
				# Format the input arguments.

  print(MFILE
	"function $outargs$fname($inargs)\n",
	"% $outargs$fname($inargs)\n",
	"$doc",
	"  $outargs$module_name(",
	join(", ", scalar(@function_names)-1, @inargs, @{$faa->{modifies}}),
	");\n");
				# Call the main MEX functions.
  close(MFILE);			# Done with that file.

  "}\n\n";			# Close off the opening brace.
}

#
# get_c_arg_scalar(\%function_def, $argname)
#
# Returns C statements to load the current value of the given argument
# into the C variable C<$function_def{args}{$argname}{c_var_name}>.  The
# variable is guaranteed to be either a scalar or an array with dimensions
# 1,1,1....
#
sub get_c_arg_scalar {
  my ($faa, $argname) = @_; # Name the arguments.

  my $arg = $faa->{args}{$argname}; # Access the definition of this argument.

  my $argtype = $arg->{basic_type}; # Get the basic type, without any frills.

  if (exists($typemap_get_scalar{$argtype})) { # Do we understand this type?
    return &{$typemap_get_scalar{$argtype}}($arg, $argname, $argtype); # Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return ("  if (!_get_pointer((void **)&$arg->{c_var_name}, $arg->{mat_expr}, @{[pointer_type_code($argtype)]})) /* $argtype */\n" .
	    "    mexErrMsgTxt(\"Wrong pointer type, expecting $argtype for argument $argname\");\n");
				# Treat it as an array of pointers.
  } else {			# Unrecognized type?
    die("wrap_matlab: don't understand type '$argtype' as scalar\n");
  }
}

#
# get_c_arg_ptr(\%function_def, $argname)
#
# Returns C statements to set up a pointer which points to the first
# value in the given argument.  The dimensions are guaranteed to be
# correct.  The type of the argument should be checked.  The pointer
# value should be stored in the variable
# $function_def{args}{$argname}{c_var_name}.
#
# The pointer should actually point to the array of all the values of
# the variable.  The array should have the same number of elements as
# the argument, since to vectorize the function, the wrapper function
# will simply step through this array.  If we want a float type and the
# input vector is double or int, then a temporary array must be made
# which is a copy of the double/int arrays.
#
sub get_c_arg_ptr {
  my ($faa, $argname) = @_; # Name the arguments.

  my $arg = $faa->{args}{$argname}; # Access the definition of this argument.

  my $argtype = $arg->{basic_type}; # Get the basic type, without any frills.

  if (exists($typemap_get_ptr{$argtype})) { # Do we understand this type?
    return &{$typemap_get_ptr{$argtype}}($arg, $argname, $argtype); # Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return ("  $arg->{c_var_name} =\n" .
	    "    MEM_ALLOC($argtype *, _arraylen($arg->{mat_expr}), sizeof ($argtype));\n" .
	    "  if (!_get_pointer((void **)$arg->{c_var_name}, $arg->{mat_expr}, @{[pointer_type_code($argtype)]})) /* $argtype */\n" .
	    "    mexErrMsgTxt(\"Wrong pointer type, expecting $argtype for argument $argname\");\n");
				# Treat it as an array of pointers.
  } else {			# Unrecognized type?
    die("wrap_matlab: don't understand type '$argtype' as scalar/matrix\n");
  }
}

#
# Get the name of the output file given the input files.  Arguments:
# 1) A reference to a list of input files.
#
sub get_outfile {
  die("wrap_matlab: must explicitly specify the output file name with -o fname.c\n");
}

#
# get_size(\%function_def, $argname, $n)
#
# Returns a C expression which is the size of the n'th dimension of the given
# argument.  Dimension 0 is the least-significant dimension.
#
sub get_size {
  my ($faa, $argname, $index) = @_; # Name the arguments.

  my $arg = $faa->{args}{$argname};

  "_dim($arg->{mat_expr}, $index)";
}


#
# Initialize the output file.  Arguments:
# 1) The name of the output file.
# 2) A reference to a list of input files explicitly listed.
# 3) A reference to the words passed to the C preprocessor.
# 4) A string that represents our guess as to the #includes which are
#    required.
#
sub initialize {
  my ($arg_outfile, $infiles, $cpp_cmd, $include_str, $ccode_str) = @_; # Name the arguments.

  $outfile = $arg_outfile;	# Store the output file name in a global variable.
  open(OUTFILE, "> $outfile") ||
    die("wrap_matlab: can't open output file $outfile--$!\n");

  print OUTFILE "/*
 * This file was automatically generated for MATLAB by wrap_matlab on
 * ", scalar(localtime(time)), "
 * from @{[@$infiles, @$cpp_cmd]}.
 */

$include_str
";

  if ($use_mxCalloc) {		# Use Matlab's allocator?
    print OUTFILE "#define MEM_ALLOC(type, n, size) ((type)mxCalloc(n, size))\n";
  } else {
    print OUTFILE "#ifndef __GNUC__
#include <alloca.h>
#endif
#define MEM_ALLOC(type, n, size) ((type)alloca((n)*(size)))

";
  }

  while (defined($_ = <DATA>)) { # Read the stuff at the end of this file.
    print OUTFILE $_;
  }

  print OUTFILE $ccode_str;

#
# We're also probably going to have to make a bunch of auxiliary files, too,
# so we can store documentation strings and also so we can write matlab
# stub functions to supply the additional function ID argument.
#
  if ($outfile =~ s@^(.*)/@@) {	# Was a directory specified for the output?
    $outdir = $1		# Use the same directory.
      unless defined($outdir);	# Unless an output directory was specified.
  } else {
    $outdir = '.'		# Just put them in the current directory.
      unless defined($outdir);
  }

  $outfile =~ s/\..*$//;	# Strip off the extension too.  (We just 
				# stripped off the directory path.)
  $module_name = $outfile;	# Remember the name of this module.  This will
				# be the name of the main MEX function.
  return "matlab::OUTFILE";	# Return the file handle.
}

#
# make_output_scalar(\%function_def, $argname)
#
#   Return C code to create the given output variable, which is guaranteed
#   to be a scalar.
#
sub make_output_scalar {
  my ($faa, $argname) = @_;
  my $arg = $faa->{args}{$argname}; # Access the definition of this argument.

  my $argtype = $arg->{basic_type}; # Get the basic type, without any frills.

  if (exists($typemap_output_scalar_make{$argtype})) { # Do we understand this type?
    return &{$typemap_output_scalar_make{$argtype}}($arg, $arg->{mat_expr});
				# Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return "  $arg->{mat_expr} = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);\n";
  } else {			# Unrecognized type?
    die("wrap_matlab: don't understand scalar output/modify type '$argtype'\n");
  }
}

#
#  make_output_ptr(\%function_def, $argname, $n_dimensions, @dimensions)
#
# Return C code to set up the given output variable.  $n_dimensions is a C
# expression, not a constant.  @dimensions is a list of C expressions that are
# the sizes of each dimension.  There may be more values in @dimensions than 
# are wanted.
#
sub make_output_ptr {
  my ($faa, $argname, $n_dims, @dims) = @_;

  my $arg = $faa->{args}{$argname}; # Access the argument info.
  my $type = $arg->{basic_type}; # Access the type.

  my $retstr = '';

  push(@dims, 1) while (@dims < 2); # Make the first two dimensions 1 if it's
				# a scalar, since that's the way matlab
				# represents it.

#
# Make an array containing the appropriate dimensions:
#
  $retstr .= ("  {\n" .
	      "    int _ds[] = { " . join(", ", @dims) . " };\n");

#
# Now make an argument of the appropriate type:
#
  if (exists($typemap_output_array_make{$type})) { # Do we understand this type?
    $retstr .= &{$typemap_output_array_make{$type}}($arg, $argname, $n_dims, "_ds");
  } elsif ($type =~ /\*$/) {	# Some random pointer type?
    $retstr .= ("    $arg->{mat_expr} = matwrapCreateNumericArray(max($n_dims, 2), _ds, mxDOUBLE_CLASS, mxCOMPLEX);\n" .
				# Allocate an array to store the result.
		"    $arg->{c_var_name} = MEM_ALLOC($type *, _arraylen($arg->{mat_expr}), sizeof ($type));\n");
				# Allocate a temporary array to put ptrs into.
  } else {
    die("wrap_matlab: can't handle scalar/matrix output of type '$type'\n");
  }

  $retstr .= "  }\n";		# Close off the definition.

  $retstr;
}

#
# n_dimensions(\%function_def, $argname)
#
# Returns a C expression which is the number of dimensions of the argument 
# whose name is $argname.
#
sub n_dimensions {
  my ($faa, $argname) = @_; # Name the arguments.

  my $arg = $faa->{args}{$argname};
  if ($arg->{basic_type} eq 'char *' || $arg->{basic_type} eq 'm_Object') {
    '0';			# Strings are treated as atomic objects,
				# and we don't support arrays of strings.
  } else {
    "_n_dims($arg->{mat_expr})";
  }
}

#
# Parse the argument vector for language-specific options.  Arguments:
# 1) A reference to the argument vector.
# Any arguments we use should be removed from the argument vector.
#
sub parse_argv {
  my ($ARGV) = @_;		# Access the arguments.
  my $argidx;

  for ($argidx = 0; $argidx < @$ARGV; ++$argidx) { # Look at the arguments:
    if ($ARGV->[$argidx] eq '-outdir') { # Output directory specified?
      $outdir = $ARGV->[$argidx+1]; # Use the specification.
      splice(@$ARGV, $argidx, 2); # Remove these two elements from the array.
      $argidx--;		# Back up to account for the arg we removed.
    } elsif ($ARGV->[$argidx] =~ /^-use_mxCalloc$/i) { # Use Matlab's memory manager.
      $use_mxCalloc = 1;	# Remember this.
      splice(@$ARGV, $argidx, 1); # Remove this argument from the list.
    }
  }
}

#
# Returns code to convert to and from pointer types to the languages
# internal representation, if any special code is needed.  If this
# subroutine is not called, then there are no class types and pointers
# will not need to be handled.
#
sub pointer_conversion_functions {
  ($_ = <<"###END") =~ s/^# ?//mg; # Strip the '# ' from beginning of lines.
# //
# // Local functions dealing with pointers:
# //
# // Pointers are stored as complex numbers.  The real part is the pointer
# // value, and the complex part is the type code.  Since the real part will be
# // a double precision number, the pointer is guaranteed to fit, even on
# // 64-bit machines.
# //
# // Fill an array with pointers.  Arguments:
# // 1) A pointer to the array to fill.  The array is already the correct size.
# // 2) The matlab object.
# // 3) The pointer type code.
# //
# // Returns false if some of the pointers have the wrong type.
# //
# static int
# _get_pointer(void **arr, const mxArray *mat, unsigned typecode)
# {
#   int idx;
#   int n_elements;
#   double *pr, *pi;
#
#   if (!mxIsComplex(mat))	// Make sure it's the right type.
#     return 0;
#
#   n_elements = _arraylen(mat);	// Get the number of elements to convert.
#   pr = mxGetPr(mat);		// Get a pointer to the real part.
#   pi = mxGetPi(mat);		// Get a pointer to the complex part.
#   for (idx = 0; idx < n_elements; ++idx)
#   {				// Look at each element:
#     unsigned this_typecode = (unsigned)*pi++; // Get the type code.
#     if (this_typecode == typecode) // Correct type?
#       *arr++ = *(void **)pr;	// Do the conversion.
#     else
#     {
#       *arr = __cvt_type(*(void **)pr, this_typecode, typecode);
# 				// Check inheritance.
#       if (*arr++ == 0) return 0; // Wrong type.
#     }
#
#     ++pr;			// Point to the next element.
#   }
#
#   return 1;			// No error.
# }
#
# //
# // Load a matrix with a pointer or an array of pointers.  This assumes
# // that the matrix has already been allocated to the correct size
# // (that's how we know how many elements to convert).
# //
# static void
# _put_pointer(mxArray *mat, void **ptr_array, double type_code)
# {
#   int idx;
#   int n_elements;
#   double *pr, *pi;
#  
#   n_elements = _arraylen(mat);	// Get the number of elements to convert.
#   pr = mxGetPr(mat);		// Get a pointer to the real part.
#   pi = mxGetPi(mat);		// Get a pointer to the complex part.
#
#   for (idx = 0; idx < n_elements; ++idx)
#   {
#     *pi++ = type_code;		// Remember the type code.
#     *pr = 0;			// Wipe out any extra bits the address
# 				// doesn't cover, in case we're on a
# 				// 32 bit machine.
#     *(void **)pr = *ptr_array++; // Remember the address.
#     ++pr;			// Advance the pointer.
#   }
# }
###END
  $_;
}

#
# put_val_scalar(\%function_def, $argname)
#
# Returns C code to take the value from the C variable whose name is given
# by $function_def{args}{$argname}{c_var_name} and store it back in the
# scripting language scalar variable.
#
sub put_val_scalar {
  my ($faa, $argname) = @_;
  my $arg = $faa->{args}{$argname}; # Access the definition of this argument.

  my $argtype = $arg->{basic_type}; # Get the basic type, without any frills.

  if (exists($typemap_put_scalar{$argtype})) { # Do we understand this type?
    return &{$typemap_put_scalar{$argtype}}($arg, $arg->{mat_expr});
				# Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return "  _put_pointer($arg->{mat_expr}, (void **)&$arg->{c_var_name}, @{[pointer_type_code($argtype)]}); /* $argtype */\n";
  } else {			# Unrecognized type?
    die("wrap_matlab: don't understand scalar output/modify type '$argtype'\n");
  }
}

#
# put_val_ptr(\%function_def, $argname)
#
# Returns C code to take the value from the C array whose name is given by
# $function_def{args}{$argname}{c_var_name} and store it back in the scripting
# language array at the specified index.  The pointer
# $function_def{args}{$argname}{c_var_name} was set up by either get_c_arg or
# make_output, depending on whether this is an input/modify or an output
# variable.
#
sub put_val_ptr {
  my ($faa, $argname) = @_;

  my $arg = $faa->{args}{$argname}; # Access the definition of this argument.

  my $argtype = $arg->{basic_type}; # Get the basic type, without any frills.
  my $argdim = $arg->{dimension}; # Point to the dimensions.

  if (exists($typemap_put_ptr{$argtype})) { # Do we understand this type?
    return &{$typemap_put_ptr{$argtype}}($arg);	# Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return "  _put_pointer($arg->{mat_expr}, (void **)$arg->{c_var_name}, @{[pointer_type_code($argtype)]}); /* $argtype */\n";
  } else {			# Unrecognized type?
    die("wrap_matlab: don't understand scalar/vector output/modify type '$argtype'\n");
  }
}


###############################################################################
#
# Various typemaps:
#

#
# Get the value of an argument which is guaranteed to be a scalar, and put
# it into a C variable.
#
# All functions in this hash are called with three arguments:
# 1) The $function_def{args}{$argname} associative array.
#    The C variable the value is stored in is $args->{c_var_name}.
# 2) The name of the argument.
# 3) The type of the argument.
#

# DSB: m_Object getter
$typemap_get_scalar{'m_Object'} = sub {
    my ($arg, $argname, $argtype) = @_; # Name the arguments.
    ("  $arg->{c_var_name} = (mxArray*) $arg->{mat_expr};\n");
};

$typemap_get_scalar{'double'} = $typemap_get_scalar{'float'} = 
  $typemap_get_scalar{'int'} = $typemap_get_scalar{'unsigned'} = 
  $typemap_get_scalar{'unsigned short'} = $typemap_get_scalar{'short'} = sub {
    my ($arg, $argname, $argtype) = @_; # Name the arguments.
    ("  if (!_get_numeric($arg->{mat_expr}, &$arg->{c_var_name}))\n" .
     "    mexErrMsgTxt(\"Expecting numeric scalar for argument $argname\");\n");
};

$typemap_get_scalar{'char'} = sub {
  my ($arg, $argname) = @_;	# Name the arguments.
  ("  if (!mxIsChar($arg->{mat_expr}))\n" .
   "    mexErrMsgTxt(\"Expecting string for argument $argname\");\n" .
   "  {\n" .
   "    char *_junk_str = MEM_ALLOC(char *, mxGetN($arg->{mat_expr}) + 1, sizeof (char));\n" .
   "    mxGetString($arg->{mat_expr}, _junk_str, mxGetN($arg->{mat_expr})+1);\n" .
   "    $arg->{c_var_name} = *_junk_str;\n" . # Just take the first char.
   "  }\n");
};

$typemap_get_scalar{'char *'} = sub {
  my ($arg, $argname) = @_;	# Name the arguments.
  ("  if (!mxIsChar($arg->{mat_expr}))\n" .
   "    mexErrMsgTxt(\"Expecting string for argument $argname\");\n" .
   "  $arg->{c_var_name} = MEM_ALLOC(char *, mxGetN($arg->{mat_expr}) + 1, sizeof (char));\n" .
   "  mxGetString($arg->{mat_expr}, $arg->{c_var_name}, mxGetN($arg->{mat_expr})+1);\n");
};

$typemap_get_scalar{'complex < double >'} =
  $typemap_get_scalar{'complex < float >'} = sub {
    my ($arg, $argname, $argtype) = @_;	# Name the arguments.
    ("  if (_get_complex($arg->{mat_expr}, &$arg->{c_var_name}) == 0)\n" .
     "    mexErrMsgTxt(\"Expecting numeric scalar for complex argument $argname\");\n");
  };

#
# Get the value of a potentially vector argument and put it into a C variable.
# The dimensions are guaranteed to be proper.  Type checking should be done.
# All functions in this hash are called with three arguments:
# 1) The $function_def{args}{$argname} associative array.
#    The C variable the value is stored in is $args->{c_var_name}.
# 2) The name of the argument.
# 3) The type of the argument.
#

$typemap_get_ptr{'double'} =	# For double precision, if the input array is
				# of Real type, we can just return a pointer
				# into it.  If it's an integer array, we have
				# to copy it into a temporary array of double.
  sub {
    my ($arg, $argname) = @_;	# Name the arguments.
    ("  if (mxIsDouble($arg->{mat_expr}))\n" .
     "    $arg->{c_var_name} = mxGetPr($arg->{mat_expr});\n" .
     "  else\n" .
     "  {\n" .
     "    $arg->{c_var_name} = MEM_ALLOC(double *, _arraylen($arg->{mat_expr}), sizeof (double));\n" .
     "    if (!_get_numeric($arg->{mat_expr}, $arg->{c_var_name}))\n" .
     "      mexErrMsgTxt(\"Expecting numeric matrix for argument $argname\");\n" .
     "  }\n");
  };

$typemap_get_ptr{'float'} =	# Float types always require a temporary array
				# since matlab doesn't have any internal
				# representation that corresponds to a float.
  sub {
    my ($arg, $argname) = @_;	# Name the arguments.
    ("  $arg->{c_var_name} = MEM_ALLOC(float *, _arraylen($arg->{mat_expr}), sizeof (float));\n" .
     "  if (!_get_numeric($arg->{mat_expr}, $arg->{c_var_name}))\n" .
     "    mexErrMsgTxt(\"Expecting numeric matrix for argument $argname\");\n");
  };

$typemap_get_ptr{'int'} = $typemap_get_ptr{'unsigned'} =
  $typemap_get_ptr{'short'} = $typemap_get_ptr{'unsigned short'} =
				# Must be an integer or an integer array.
  sub {
    my ($arg, $argname, $argtype) = @_;	# Name the arguments.
    my $ilen = ($argtype =~ /\bshort$/ ? '16' : '32');
				# Which Matlab subroutines to call to check
				# for 16 or 32 bit numbers.
    ("  if (mxIsInt$ilen($arg->{mat_expr}) || mxIsUint$ilen($arg->{mat_expr}))\n" .
				# Is it already an integer?
     "    $arg->{c_var_name} = ($argtype *)mxGetData($arg->{mat_expr});\n".
     "  else\n" .
     "  {\n" .
     "    $arg->{c_var_name} = MEM_ALLOC($argtype *, _arraylen($arg->{mat_expr}), sizeof ($argtype));\n" .
     "    if (!_get_numeric($arg->{mat_expr}, $arg->{c_var_name}))\n" .
     "      mexErrMsgTxt(\"Expecting numeric matrix for argument $argname\");\n" .
     "  }\n");
  };

$typemap_get_ptr{'complex < float >'} =	# Complex numbers require temporary 
  $typemap_get_ptr{'complex < double >'} = sub { # arrays because matlab
				# stores the real part in a different location
				# from the complex part.
    my ($arg, $argname, $argtype) = @_;	# Name the arguments.

    ("  $arg->{c_var_name} = MEM_ALLOC($argtype *, _arraylen($arg->{mat_expr}), sizeof ($argtype));\n" .
     "  if (!_get_complex($arg->{mat_expr}, $arg->{c_var_name}))\n" .
     "    mexErrMsgTxt(\"Expecting complex matrix for argument $argname\");\n");
};

$typemap_get_ptr{'char *'} = sub {
  die("wrap_matlab: arrays of strings are not yet supported\n");
};	# Arrays of strings are not yet supported.

$typemap_get_ptr{'char'} = $typemap_get_ptr{'unsigned char'} = sub {
  my ($arg, $argname, $argtype) = @_;	# Name the arguments.
  ("  $arg->{c_var_name} = MEM_ALLOC($argtype *, _arraylen($arg->{mat_expr})+1, sizeof ($argtype));\n" .
				# Allocate a temporary array for the string.
   "  if (!mxIsChar($arg->{mat_expr}))\n" .
   "    mexErrMsgTxt(\"Expecting string for argument $argname\");\n" .
   "  mxGetString($arg->{mat_expr}, $arg->{c_var_name}, mxGetN($arg->{mat_expr})+1);\n");
};

#
# Create an output argument vector of a given size.  Arguments:
# 1) The argument description associative array.
# 2) The argument name.
# 3) The number of dimensions.  Since the number of dimensions is not
#    necessarily known in advance, this is a C expression which gives the
#    number of dimensions.
# 4) The name of the integer array containing the structure.
#
$typemap_output_array_make{'double'} = 
  $typemap_output_array_make{'float'} = # All of these types are returned as
  $typemap_output_array_make{'int'} = # double precision numbers.
  $typemap_output_array_make{'unsigned'} =
  $typemap_output_array_make{'short'} =
  $typemap_output_array_make{'unsigned short'} = sub {
    my ($arg, $argname, $n_dims, $c_dim_arr) = @_; # Name the arguments.

    my $retstr = "    $arg->{mat_expr} = matwrapCreateNumericArray(max($n_dims, 2), $c_dim_arr, mxDOUBLE_CLASS, mxREAL);\n";
				# Create space to hold the result.
    if ($arg->{basic_type} eq 'double') { # No temporary array required?
      $retstr .= "    $arg->{c_var_name} = mxGetPr($arg->{mat_expr});\n";
				# Access the array of doubles directly.
    } else {
      $retstr .= "    $arg->{c_var_name} = MEM_ALLOC($arg->{basic_type} *, _arraylen($arg->{mat_expr}), sizeof ($arg->{basic_type}));\n";
    }

    $retstr;
  };

$typemap_output_array_make{'char'} =
  $typemap_output_array_make{'unsigned char'} = sub {
    my ($arg, $argname, $n_dims, $c_dim_arr) = @_; # Name the arguments.
    
    ("    $arg->{mat_expr} = mxCreateCharArray(max($n_dims, 2), $c_dim_arr);\n" .
				# Make the array.
     "    $arg->{c_var_name} = MEM_ALLOC($arg->{basic_type} *, _arraylen($arg->{mat_expr}), sizeof ($arg->{basic_type}));\n");
				# Allocate temporary space for the string.
				# Matlab stores strings as uint16 characters,
				# which is definitely not what we want to use.
  };

$typemap_output_array_make{'complex < float >'} =
  $typemap_output_array_make{'complex < double >'} = sub {
    my ($arg, $argname, $n_dims, $c_dim_arr) = @_; # Name the arguments.
    
    ("    $arg->{mat_expr} = matwrapCreateNumericArray(max($n_dims, 2), $c_dim_arr, mxDOUBLE_CLASS, mxCOMPLEX);\n" .
				# Make the array.
     "    $arg->{c_var_name} = MEM_ALLOC($arg->{basic_type} *, _arraylen($arg->{mat_expr}), sizeof ($arg->{basic_type}));\n");
				# Allocate space for a temporary array.
  };

#
# Create a scalar value for output.  Arguments:
# 1) The associative array describing the argument.
# 2) A C expression for the argument.

# DSB: m_Object -- nothing to do
$typemap_output_scalar_make{'m_Object'} = sub { ''; };

$typemap_output_scalar_make{'double'} = $typemap_output_scalar_make{'float'} =
  $typemap_output_scalar_make{'int'} =
  $typemap_output_scalar_make{'unsigned'} = sub {
    "  $_[1] = mxCreateDoubleMatrix(1, 1, mxREAL);\n"; # Create the matrix.
  };

$typemap_output_scalar_make{'complex < double >'} =
  $typemap_output_scalar_make{'complex < float >'} = sub {
    "  $_[1] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);\n";
  };

$typemap_output_scalar_make{'char'} = sub {
  if (!defined($output_scalar_dimensions)) { # Did we output the scalar
				# dimension constants?
    print OUTFILE "
static int scalar_dimensions[] = {1,1}; // A value to pass to give the 
				// dimensions of a scalar.

";
				# Output with print so it goes before the
				# function we are wrapping.
    $output_scalar_dimensions = 1; # Don't output it again.
  }


  "  $_[1] = mxCreateCharArray(1, scalar_dimensions);\n";
};

$typemap_output_scalar_make{'char *'} = sub { ''; };
				# Nothing to do for string scalars.

#
# Set the value of an output/modify argument which is guaranteed to be a
# scalar.  All functions in this hash are called with three arguments:
# 1) The $function_def{args}{$argname} associative array.
# 2) A C expression for the argument.
#

# DSB: m_Object
$typemap_put_scalar{'m_Object'} = sub {
    my ($arg, $matexpr) = @_;	# Name the arguments.
    "  $matexpr = $arg->{c_var_name};\n";
  };

$typemap_put_scalar{'double'} = $typemap_put_scalar{'float'} =
  $typemap_put_scalar{'int'} = $typemap_put_scalar{'unsigned'} =
  $typemap_put_scalar{'short'} = $typemap_put_scalar{'unsigned short'} = sub {
    my ($arg, $matexpr) = @_;	# Name the arguments.
    "  *mxGetPr($matexpr) = (double)$arg->{c_var_name};\n";
  };

$typemap_put_scalar{'complex < double >'} =
  $typemap_put_scalar{'complex < float >'} = sub {
    my ($arg, $mat_expr) = @_;	# Name the arguments.
    ("  *mxGetPr($mat_expr) = (double)$arg->{c_var_name}.real();\n" .
     "  *mxGetPi($mat_expr) = (double)$arg->{c_var_name}.imag();\n");
  };

$typemap_put_scalar{'char'} = sub {
  my ($arg, $mat_expr) = @_;	# Name the arguments.
  "  *(char *)mxGetPr($mat_expr) = $arg->{c_var_name};\n";
};

$typemap_put_scalar{'char *'} = sub {
  my ($arg, $mat_expr) = @_;	# Name the arguments.
  "  $mat_expr = mxCreateString($arg->{c_var_name});\n";
};

#
# Set the value of an output/modify argument which may be an array.
# The type and dimensions are guaranteed to be proper.
# All functions in this hash are called with three arguments:
# 1) The $function_def{args}{$argname} associative array.
#
$typemap_put_ptr{'double'} = sub { "" }; # Nothing to do for a double, since
				# we already manipulated the data in the
				# matlab object.

$typemap_put_ptr{'float'} =	# These are more complicated because we need
  $typemap_put_ptr{'int'} =	# the intermediate array.
  $typemap_put_ptr{'unsigned'} =
  $typemap_put_ptr{'short'} =
  $typemap_put_ptr{'unsigned short'} = sub {
    my ($arg) = @_; # Name the arguments.

    ("  {\n" .
     "    int _fidx;\n" .
     "    double *d_ptr = mxGetPr($arg->{mat_expr});\n" . # Point to output arr.
     "    for (_fidx = 0; _fidx < _arraylen($arg->{mat_expr}); ++_fidx)\n" .
     "      *d_ptr++ = (double)*$arg->{c_var_name}++;\n" . # Copy each element.
     "  }\n");
  };

$typemap_put_ptr{'complex < double >'} = $typemap_put_ptr{'complex < float >'}
= sub {
  my ($arg) = @_; # Name the arguments.

  ("  {\n" .
   "    int _fidx;\n" .
   "    double *r_ptr = mxGetPr($arg->{mat_expr});\n" . # Point to output arr.
   "    double *i_ptr = mxGetPi($arg->{mat_expr});\n" . # Point to output arr.
   "    for (_fidx = 0; _fidx < _arraylen($arg->{mat_expr}); ++_fidx)\n" .
   "    {\n" .
   "      *r_ptr++ = (double)$arg->{c_var_name}" . "->real();\n" .
   "      *i_ptr++ = (double)$arg->{c_var_name}" . "->imag();\n" .
   "      ++$arg->{c_var_name};\n" . # Go to next position.
   "    }\n" .
   "  }\n");
};
  
$typemap_put_ptr{'char'} = $typemap_put_ptr{'unsigned char'} = sub {
  my ($arg) = @_; # Name the arguments.
  ("  {\n" .
   "    int _cidx;\n" .
   "    mxChar *mc_ptr = (mxChar *)mxGetPr($arg->{mat_expr});\n" .
   "    for (_cidx = 0; _cidx < _arraylen($arg->{mat_expr}); ++_cidx)\n" .
   "      *mc_ptr++ = (mxChar)*$arg->{c_var_name}++;\n" .
   "  }\n");
};

__DATA__
#ifdef __cplusplus
extern "C" {
#endif
#include "matrix.h"
#include "mex.h"
#ifdef __cplusplus
}
#endif
#include <complex>
#define complex std::complex

//
// Define m_Object for pass-through of MATLAB types:
//
typedef mxArray* m_Object;

//
// Reserve space for first_arg (used as a hack for getting mwSize
// if needed -- would just use mwSize if it didn't break compatibility
// with pre-7.3 MATLAB).
// 
const mxArray* first_arg;

//
// Wrapper around mxCreateNumericArray for going from int to mwSize,
// if necessary.
// 
#define matwrapCreateNumericArray(ndim, dims, classid, complex_flag) \
   matwrapCreateNumericArray1(ndim, dims, classid, complex_flag, \
                              mxGetM(first_arg))

template <class S>
inline mxArray* 
matwrapCreateNumericArray1(int ndim, int* dims, mxClassID classid,
                           mxComplexity complex_flag, S dummy)
{
    S dims1[20];
    for (int i = 0; i < ndim; ++i)
        dims1[i] = dims[i];
    return mxCreateNumericArray(ndim, dims1, classid, complex_flag);
}

//
// Return the maximum of the arguments:
//
template <class FLOAT>
inline FLOAT
max(FLOAT arg1, FLOAT arg2)
{
  return arg1 >= arg2 ? arg1 : arg2;
}

/*
 * Return the number of dimensions of an array.
 * This is not well-defined in matlab; we just return the highest dimension
 * which is not 1.
 */
static int
_n_dims(const mxArray *mat)
{
  int n_dim = mxGetNumberOfDimensions(mat); /* Get what matlab thinks. */
  if (n_dim == 2 && mxGetN(mat) == 1) /* 2D array with # of columns = 1? */
  {
    n_dim = 1;			/* It's really a 1D array. */
    if (mxGetM(mat) == 1)	/* Also only one row? */
      n_dim = 0;		/* It's really a scalar. */
  }

  return n_dim;
}

/*
 * Function to get the total number of elements in an array.
 * Arguments:
 * 1) The matlab object.
 *
 * If the object is a scalar, the array length is 1.
 */
static int
_arraylen(const mxArray *mat) {
  int len = 1;
  int idx;
  for (idx = 0; idx < mxGetNumberOfDimensions(mat); ++idx)
    len *= mxGetDimensions(mat)[idx];		// Multiply all the dimensions.

  return len;
}

//
// Get the n'th dimension of a Matlab object.
// Arguments:
// 1) The matlab matrix.
// 2) Which dimension.  0 for the most quickly varying one, etc.
//
static inline int
_dim(const mxArray *mat, int idxno) {
  if (idxno >= mxGetNumberOfDimensions(mat))	// Illegal # of dims?
    return 1;
  else
    return mxGetDimensions(mat)[idxno];	// Return the dimension.
}

//
// Load up an array (which may consist only of a single element) with
// the values from a matlab object.  Arguments:
// 1) The matlab object.
// 2) Where to put the values.  The array should be allocated sufficiently
//    large so that all of the values in the matlab object can be converted.
//
// These functions return 0 on error, or 1 if the operation succeded.
//
template <class FLOAT>
static int
_get_numeric(const mxArray *mat, FLOAT *arr) // float or double.
{
  int idx;
  int n_elements = _arraylen(mat); // Get the length.

  if (mxIsDouble(mat))		// Ordinary double precision matrix?
  {
    double *pr = mxGetPr(mat);	// Point to vector of double precision values.

    for (idx = 0; idx < n_elements; ++idx) // Loop through all elements.
      *arr++ = (FLOAT)*pr++;	// Convert and copy this element.
  }
  else if (mxIsInt32(mat) || mxIsUint32(mat)) // 32 bit integer?
  {
    INT32_T *pr = (INT32_T *)mxGetData(mat); // Point to the data.
    for (idx = 0; idx < n_elements; ++idx) // Loop through all elements.
      *arr++ = (FLOAT)*pr++;	// Convert and copy this element.
  }
  else if (mxIsInt16(mat) || mxIsUint16(mat)) // 16 bit integer?
  {
    INT16_T *pr = (INT16_T *)mxGetData(mat); // Point to the data.
    for (idx = 0; idx < n_elements; ++idx) // Loop through all elements.
      *arr++ = (FLOAT)*pr++;	// Convert and copy this element.
  }
  else if (mxIsInt8(mat) || mxIsUint8(mat)) // 8 bit integer?
  {
    INT8_T *pr = (INT8_T *)mxGetData(mat); // Point to the data.
    for (idx = 0; idx < n_elements; ++idx) // Loop through all elements.
      *arr++ = (FLOAT)*pr++;	// Convert and copy this element.
  }
  else				// Don't know what this is.
    return 0;			// Type conversion error.

  return 1;			// No error.
}

//
// A version of the above function for complex numbers:
//
template <class FLOAT>
static int
_get_complex(const mxArray *mat, complex<FLOAT> *arr)
{
  int idx;
  int n_elements = _arraylen(mat); // Get the length.

  if (mxIsDouble(mat))		// Ordinary double precision matrix?
  {
    double *pr = mxGetPr(mat);	// Point to vector of double precision values.
    double *pi = mxGetPi(mat);	// Get complex part, if any.

    if (pi != 0)		// Is there a complex part?
    {
      for (idx = 0; idx < n_elements; ++idx) // Loop through all elements.
        *arr++ = complex<FLOAT>((FLOAT)*pr++, (FLOAT)*pi++);
				// Convert and copy this element.
    }
    else			// Purely real matrix:
    {
      for (idx = 0; idx < n_elements; ++idx) // Loop through all elements.
        *arr++ = complex<FLOAT>((FLOAT)*pr++, 0);
				// Convert and copy this element, setting
				// the imaginary part to 0.
    }

    return 1;			// No error.
  }
  else				// Currently we don't support the int32, int16,
				// or int8 types.
    return 0;			// Type conversion error.
}
