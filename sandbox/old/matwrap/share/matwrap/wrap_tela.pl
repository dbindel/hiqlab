#
# Language module for the wrapper generator for tela.  Most of these
# functions are called directly from the wrapper generator.  Internally used
# functions are prefixed with "tela_".
#
# Copyright (c) 1997 Gary R. Holt.  This is distributed under the terms of the 
# perl artistic license (http://language.perl.com/misc/Artistic.html).
#
package tela;			# Everything should be in this package.

($progname = $0) =~ s@.*/@@;	# Get the program name for error messages.

$max_dimensions = 4;		# Maximum number of dimensions in tensors.

*pointer_type_code = *main::pointer_type_code;
				# Copy the definition of pointer_type_code
				# into this module.

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
  $argname;			# Just pass by name.
}

#
# arg_declare("arg_name_in_arglist")
# 
# This returns a C/C++ declaration appropriate for the argument passed
# using arg_pass.  For example, in the MATLAB module this function returns
# "mxArray *arg_name_in_arglist".
#
sub arg_declare {
  "const Tobject &$_[0]";
}

#
# declare_const("constant name", "class name", "<type>", "doc str")
#
# Output routines to make a given constant value accessible from the 
# interpreter.
# If "class name" is blank, this is a global constant.
#
sub declare_const {
  my ($const_name, $class, $type, $doc) = @_; # Name the arguments.
#
# We wrap constants currently by providing a simple get function.
#
  "[retval] = get_$const_name()
/* get_$const_name() returns the value of $const_name.
 */
{
  retval = $const_name;
  return 0;
}
";
				# That's all there is to it.  Use tela's
				# overloaded operator= to do the type
				# conversion.
}

#
# error_dimension(\%function_def, $argname)
#
# A C statement (including the final semicolon, if not surrounded by braces)
# which indicates that an error has occured because the dimension of argument
# $argname was wrong.
#
sub error_dimension {
  "return -2;";			# Currently we don't have any way to specify
				# which variable was in error.
}

#
# finish()
#
# Called after all functions have been wrapped, to close the output file and do
# whatever other cleanup is necessary.
#
sub finish {
  close(OUTFILE);		# Done with this file.
}


#
# Begin the definition of a function.  Arguments:
# 1) The %function_def array.
#
sub function_start {
  my $faa = $_[0];		# Access the argument.

#
# We don't support vectorization of char or char * arguments, so mark
# these as scalars:
#
  foreach (@{$faa->{inputs}}, @{$faa->{modifies}}, @{$faa->{outputs}}) {
    my $arg = $faa->{args}{$_};
    if ($arg->{basic_type} eq 'char' || $arg->{basic_type} eq 'char *') {
      $arg->{vectorize} = 0;	# This cannot be vectorized.
      $arg->{source} eq 'output' and $faa->{vectorize} = 0;
				# Function cannot be vectorized at all if it
				# takes char or char * output.
      @{$arg->{dimension}} and
	die("wrap_tela: vectors of char or char * are not supported\n");
    }
  }

  my $funcheader =		# Format the function header for the .ct file.
    sprintf("[%s] = %s(%s)",
	    join(", ", @{$faa->{outputs}}, @{$faa->{modifies}}),
				# Format the output/modify arguments.
	    $faa->{script_name} || ($faa->{class} ? $faa->{class} . "_" : "") . $faa->{name},
				# Form the name of the tela function.
	    join(", ", @{$faa->{inputs}}));
				# Format the input arguments.
  "$funcheader
/* $funcheader
   Error codes:
   -1: Wrong type of input/modify argument
   -2: Mismatched dimensions
 */\n{\n";			# Output the header.
}

#
# function_end(\%function_def)
#
# Return a string which finishes off the definition of a wrapper.
#
sub function_end {
  "  return 0;\n}\n\n";		# Just close off the opening brace.
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
  my ($faa, $argname) = @_;	# Name the arguments.
  my $arg = $faa->{args}{$argname}; # Access the definition of this argument.

  my $argtype = $arg->{basic_type}; # Get the basic type, without any frills.

  if (exists($typemap_get_scalar{$argtype})) { # Do we understand this type?
    return &{$typemap_get_scalar{$argtype}}($arg, $argname, $argtype); # Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return ("  if (!_get_pointer((void **)&$arg->{c_var_name}, $argname, @{[pointer_type_code($argtype)]})) /* $argtype */\n" .
	    "    return -1;\n");
				# Treat it as an array of pointers.
  } else {			# Unrecognized type?
    die("wrap_tela: don't understand type '$argtype' as scalar\n");
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
	    "    ($argtype *)alloca(sizeof ($argtype)*_arraylen($argname));\n" .
	    "  if (!_get_pointer((void **)$arg->{c_var_name}, $argname, @{[pointer_type_code($argtype)]})) /* $argtype */\n" .
	    "    return -1;\n");
				# Treat it as an array of pointers.
  } else {			# Unrecognized type?
    die("wrap_tela: don't understand type '$argtype' as scalar/matrix\n");
  }
}

#
# Get the name of the output file given the input files.  Arguments:
# 1) A reference to a list of input files.
#
sub get_outfile {
  my ($inarr) = @_;		# Access the arguments.

  my $outfile = $inarr->[-1];	# Take the name from the last file.
  $outfile =~ s/\.[^.]+$//;	# Strip off the extension.

  $outfile . ".ct";		# Add tela's extension.
}

#
# get_size(\%function_def, $argname, $n)
#
# Returns a C expression which is the size of the n'th dimension of the given
# argument.  Dimension 0 is the least-significant dimension.
#
sub get_size {
  my ($faa, $argname, $index) = @_; # Name the arguments.

  if ($faa->{args}{$argname}{basic_type} eq 'char *') {
				# Tela doesn't support vectors of char *, so
    return 1;			# these can't be vectorized.
  } else {
    "_dim($argname, $index)";
  }
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
  my ($outfile, $infiles, $cpp_cmd, $include_str) = @_;	# Name the arguments.

  open(OUTFILE, "> $outfile") ||
    die("$progname: can't open output file $outfile--$!\n");

  print OUTFILE "/*
 * This file was automatically generated for Tela by wrap_tela on
 * ", scalar(localtime(time)), "
 * from @{[@$infiles, @$cpp_cmd]}.
 */

$include_str
";				# Output a header.

  while (defined($_ = <DATA>)) { # Read the rest of the static functions.
    print OUTFILE $_;
  }
				
  return "tela::OUTFILE";	# Return the file handle.
}

#
# make_output_scalar(\%function_def, $argname)
#
#   Return C code to create the given output variable, which is guaranteed
#   to be a scalar.
#
sub make_output_scalar {
				# There's nothing to do at this point.  We
				# use Tela's overloaded operator = to do the
				# assignemnt when put_val is called.
  '';
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
  return '' if $type eq 'char *'; # We can't handle arrays of 
				# char *, and scalars are already handled
				# properly.

  my $retstr = '';

  my $min_dim = @{$arg->{dimension}}; # Get the minimum dimension.
  my $max_dim = 4;
  if (!$arg->{vectorize}) {	# Not being vectorized?
    $max_dim = $min_dim;
  }

  if ($n_dims =~ /^\d+$/) {	# Is the dimension known exactly?
    $min_dim = $max_dim = $n_dims; # Use it.
  }

  if ($max_dim > 0) {		# Can this be a vector?
    $retstr .= "  if ($n_dims != 0)\n" # First handle the array case.
      unless ($min_dim > 0);	# Don't bother with this line if scalars aren't
				# allowed.
#
# First we have to make a TDimPack object containing the appropriate
# dimensions.  This is a little tricky because the dimension list we have
# is specified in column-major order, whereas we want it in row major order.
#
    $retstr .= ("  {\n" .
		"    Tint _ds[] = { " . join(", ", reverse @dims) . " };\n" .
		"    TDimPack _dims(&_ds[$#dims+1-$n_dims], $n_dims);\n");

#
# Now make an argument of the appropriate type:
#
    if (exists($typemap_output_array_make{$type})) { # Do we understand this type?
      $retstr .= &{$typemap_output_array_make{$type}}($argname,
						      $arg->{c_var_name},
						      "_dims", $type);
    } elsif ($type =~ /\*$/) {	# Some random pointer type?
      $retstr .= ("    $argname.creserv(_dims);\n" . # Allocate an array of
				# complex number to store the result.
		  "    $arg->{c_var_name} = ($type *)alloca(_dims.length() * sizeof ($type));\n");
				# Allocate a temporary array to put ptrs into.
    } else {
      die("wrap_tela: can't handle vector output of type '$type'\n");
    }

    $retstr .= "  }\n";		# Close off the definition.
  }
#
# That's it for the vectors.  Now put in code for the scalar case.
#
  $retstr .= "  else\n  {\n" unless ($min_dim > 0 || $max_dim == 0);

  if ($min_dim == 0) {		# Are scalars allowed?
    if ($typemap_output_scalar_make{$type}) { # Is there something to do for this type?
      $retstr .= &{$typemap_output_scalar_make{$type}}($argname, $arg->{c_var_name}, $type);
    } elsif ($type =~ /\*$/) {	# Some random pointer type?
      $retstr .= ("    $arg->{c_var_name} = ($type *)alloca(sizeof ($type));\n" .
				# Leave space for it.
		  "    $argname.SetToVoid();\n"); # Clear any previous array value.
				# This is a flag that indicates that we want
				# a scalar object for output.
    } else {
      die("wrap_tela: can't handle scalar output type '$type'\n");
    }

  }
  $retstr .= "  }\n"		# Terminate the else clause.
    unless $min_dim > 0 || $max_dim == 0;

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

  $faa->{args}{$argname}{basic_type} eq 'char *' and return 0;
				# Tela doesn't support vectors of strings.
  "_n_dims($argname)";
}


#
# Parse the argument vector for language-specific options.  Arguments:
# 1) A reference to the argument vector.
# Any arguments we use should be removed from the argument vector.
#
sub parse_argv {
				# We don't do anything now.
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
# // Pointers are stored as double precision numbers (actually, as the real part
# // of a complex number, with the imaginary part being the type).  This
# // should work even on 64 bit machines, since a double is at least 64 bits
# // on all machines that I know of.
# //
# union d_to_p {
#   double dval;			// Double precision value.
#   void *pval;			// Pointer value.
# };

# //
# // Local functions dealing with pointers:
# //
# // Pointers are stored as complex numbers.  The real part is the pointer value,
# // and the complex part is the type code.  Since the real part will be a double
# // precision number, the pointer is guaranteed to fit, even on 64-bit machines.
# //
# // Fill an array with pointers.  Arguments:
# // 1) A pointer to the array to fill.  The array is already the correct size.
# // 2) The name of the source Tela object.
# // 3) The pointer type code.
# //
# // Returns false if some of the pointers have the wrong type.
# //
# static int
# _get_pointer(void **arr, const Tobject &obj, unsigned typecode) {
#   switch (obj.kind()) {
#   case Kcomplex:		// Scalar complex value?
#   {
#     d_to_p c_real;
#     c_real.dval = obj.ComplexValue().real(); // Get the pointer value.
#     unsigned this_typecode = unsigned(obj.ComplexValue().imag());
#     if (this_typecode == typecode)
#       *arr = c_real.pval;	// Correct type.  Just store the value.
#     else			// Type code does not match.  Check for
#       *arr = __cvt_type(c_real.pval, this_typecode, typecode); // inheritance.
#     return *arr != 0;		// Return 1 if it was valid.
#   }
#   case KComplexArray:
#   {
#     for (Tcomplex *c_ptr = obj.ComplexPtr();	// Look at each array element.
#          c_ptr < &obj.ComplexPtr()[obj.length()]; ++c_ptr)
#     {
#       d_to_p c_real;
#       c_real.dval = c_ptr->real(); // Get the pointer value.

#       unsigned this_typecode = unsigned(c_ptr->imag()); // Get type code.
#       if (this_typecode == typecode)	// Correct type?
# 	*arr++ = c_real.pval;	// Store it.
#       else
#       {				// Check for inheritance:
# 	*arr = __cvt_type(c_real.pval, this_typecode, typecode);
#         if (*arr++ == 0) return 0; // Not a base class.
#       }
#     }
#     break;
#   }
#   default:			// Some other type?
#     return 0;
#   }

#   return 1;			// No error.
# }

# //
# // Load a Tela object with a pointer or an array of pointers:
# //
# static void
# _put_pointer(Tobject &obj, void **ptr_array, double type_code) {
#   if (_n_dims(obj) == 0)	// Is this a scalar output?
#   {
#     d_to_p c_real;
#     c_real.dval = 0;		// Zero out the unused part of the field.
#     c_real.pval = *ptr_array;	// Get the value.
#     obj = Tcomplex(c_real.dval, type_code); // Set the object.
#   }
#   else				// It's a vector output.  In this case, the
#   {				// object's dimensions are already set.
#     for (int idx = 0; idx < obj.length(); ++idx)
#     {
#       d_to_p c_real;
#       c_real.dval = 0;		// Zero out the unused part of the field.
#       c_real.pval = *ptr_array++; // Get the value.
#       obj.ComplexPtr()[idx] = Tcomplex(c_real.dval, type_code);
#     }
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
    return &{$typemap_put_scalar{$argtype}}($arg, $argname, $arg->{c_var_name});
				# Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return "    _put_pointer($argname, (void **)&$arg->{c_var_name}, @{[pointer_type_code($argtype)]}); /* $argtype */\n";
  } else {			# Unrecognized type?
    die("wrap_tela: don't understand scalar output/modify type '$argtype'\n");
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
    return &{$typemap_put_ptr{$argtype}}($arg, $argname, $arg->{c_var_name});
				# Do the conversion.
  } elsif ($argtype =~ /\*$/) { # Is this a pointer class we don't understand?
    return "  _put_pointer($argname, (void **)$arg->{c_var_name}, @{[pointer_type_code($argtype)]}); /* $argtype */\n";
  } else {			# Unrecognized type?
    die("wrap_tela: don't understand scalar/vector output/modify type '$argtype'\n");
  }
}

###############################################################################
#
# Typemap stuff.  These are subroutines that do various things depending on the
# type.  You can add freely to this to define new types.
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
$typemap_get_scalar{'double'} = $typemap_get_scalar{'float'} = sub {
  my ($arg, $argname, $argtype) = @_; # Name the arguments.
  ("  if (_get_double($argname, &$arg->{c_var_name}) == 0)\n" .
   "    return -1;\n");
};

$typemap_get_scalar{'int'} = $typemap_get_scalar{'unsigned'} = sub {
  my ($arg, $argname, $argtype) = @_; # Name the arguments.
  ("  $argtype *___p_$argname = ($argtype *)_get_int_ptr($argname);\n" .
   "  if (___p_$argname == 0)\n" . # Some error?
   "    return -1;\n" .		# It's a type error.
   "  $arg->{c_var_name} = *___p_$argname;\n");
};

$typemap_get_scalar{'char *'} = sub {
  my ($arg, $argname) = @_;	# Name the arguments.
  ("  if (!$argname.HasStringFlag()) return -1;\n" .
   "  $arg->{c_var_name} = (char *)alloca(($argname.length()+1) * sizeof (char));\n" .
   "  {\n" .
   "    int idx;\n" .
   "    for (idx = 0; idx < $argname.length(); ++idx)\n" .
   "      $arg->{c_var_name}" . "[idx] = $argname.IntPtr()[idx];\n" .
   "    $arg->{c_var_name}" . "[idx] = '\\0';\n" .
   "  }\n");
};

$typemap_get_scalar{'complex < double >'} =
  $typemap_get_scalar{'complex < float >'} = sub {
    my ($arg, $argname, $argtype) = @_;	# Name the arguments.
    ("  if (_get_complex($argname, &$arg->{c_var_name}) == 0)\n" .
     "    return -1;\n");
  };

#
# Get the value of a potentially vector argument and put it into a C variable.
# The dimensions are guaranteed to be proper.  Type checking should be done.
# -1 should be returned on a type error.
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
    ("  switch ($argname.kind())\n" .
     "  {\n" .
     "    case Kreal:\n" .
     "      $arg->{c_var_name} = &(*(Tobject *)&$argname).RealRef();\n" .
     "      break;\n" .
     "    case KRealArray:\n" .
     "      $arg->{c_var_name} = $argname.RealPtr();\n" .
     "      break;\n" .
     "    default:\n" .
     "      $arg->{c_var_name} = (double *)alloca(_arraylen($argname) * sizeof (double));\n" .
     "      if (_get_double($argname, $arg->{c_var_name}) == 0)\n" .
     "        return -1;\n" .
     "      break;\n" .
     "   }\n");
  };

$typemap_get_ptr{'float'} =	# float types always require a temporary array
				# since Tela doesn't have any internal
				# representation that corresponds to a float.
  sub {
    my ($arg, $argname) = @_;	# Name the arguments.
    ("      $arg->{c_var_name} = (float *)alloca(_arraylen($argname) * sizeof (float));\n" .
     "  if (_get_double($argname, $arg->{c_var_name}) == 0)\n" .
     "    return -1;\n");
  };

$typemap_get_ptr{'int'} = $typemap_get_ptr{'unsigned'} =
				# Must be an integer or an integer array.
  sub {
    my ($arg, $argname, $argtype) = @_;	# Name the arguments.
    ("  $arg->{c_var_name} = ($argtype *)_get_int_ptr($argname);\n" . # Get a pointer.
     "  if ($arg->{c_var_name} == 0)\n" .
     "    return -1;\n");
  };

$typemap_get_ptr{'complex < float >'} = sub {
				# For single precision complex numbers, we have
				# to convert from Tela's type.
  my ($arg, $argname, $argtype) = @_;	# Name the arguments.

  ("  $arg->{c_var_name} = ($argtype *)alloca(_arraylen($argname) * sizeof ($argtype));\n" .
   "  if (_get_complex($argname, $arg->{c_var_name}) == 0)\n" .
   "    return -1;\n");
};

$typemap_get_ptr{'complex < double >'} = sub {
				# For double precision complex numbers, we can
				# use Tela's type w/o copying.
  my ($arg, $argname, $argtype) = @_;	# Name the arguments.

  ("  switch ($argname.kind())\n" .
   "  {\n" .
   "  case Kcomplex:\n" .	# Scalar complex value?
   "    $arg->{c_var_name} = ($argtype *)&(*(Tobject *)&$argname).ComplexRef();\n" .
   "    break;\n" .
   "  case KComplexArray:\n" .	# Vector complex value?
   "    $arg->{c_var_name} = ($argtype *)$argname.ComplexPtr();\n" . # Just point to tela's object.
   "    break;\n" .
   "  default:\n" .
   "    $arg->{c_var_name} = ($argtype *)alloca(_arraylen($argname) * sizeof ($argtype));\n" .
   "    if (_get_complex($argname, $arg->{c_var_name}) == 0)\n" .
   "      return -1;\n" .
   "    break;\n" .
   "  }\n");
};

#
# Set the value of an output/modify argument which is guaranteed to be a
# scalar.  All functions in this hash are called with three arguments:
# 1) The $function_def{args}{$argname} associative array.
# 2) The name of the argument.
# 3) A C expression which is the new value of the variable.
#
$typemap_put_scalar{'double'} = $typemap_put_scalar{'float'} =
  $typemap_put_scalar{'int'} = 
  $typemap_put_scalar{'char *'} = $typemap_put_scalar{'char'} = sub {
    "  $_[1] = $_[2];\n";	# Use Tela's overloaded operator=.
  };

$typemap_put_scalar{'unsigned'} = sub {	# Can't use overloaded = here because
  "  $_[1] = int($_[2]);\n";	# tela doesn't know about unsigned integers.
};

$typemap_put_scalar{'complex < double >'} =
  $typemap_put_scalar{'complex < float >'} = sub {
    "  $_[1] = Tcomplex($_[2].real(), $_[2].imag());\n";
  };

#
# Set the value of an output/modify argument which may be an array.
# The type and dimensions are guaranteed to be proper.
# All functions in this hash are called with three arguments:
# 1) The $function_def{args}{$argname} associative array.
# 2) The name of the argument.
# 3) A C expression which is the new value of the variable.
#
$typemap_put_ptr{'double'} =	# Nothing to do for a double.
  $typemap_put_ptr{'int'} = $typemap_put_ptr{'unsigned'} = # Nor for an int.
  $typemap_put_ptr{'complex < double >'} = # Or for double precision complex.
  sub { ""; };

$typemap_put_ptr{'float'} =	# A float is more complicated because we need
  sub {				# the intermediate array.
    my ($arg, $argname, $cvar) = @_; # Name the arguments.
    my $retstr = '';
    my $could_be_scalar = !@{$arg->{dimension}}; # Non-zero if this is definitely a vector.
    my $could_be_vec = @{$arg->{dimension}} || $arg->{vectorize}; # True if this might be a vector.
    if ($could_be_vec) {	# Could this be a vector?
      $retstr .=  "  if ($argname.IsArray())\n"
	if $could_be_scalar;	# Check (unless we know for sure).
      $retstr .= ("  {\n" .
		  "    for (int _fidx = 0; _fidx < $argname.length(); ++_fidx)\n" .
		  "      $argname.RealPtr()[_fidx] = double(*$cvar++);\n" .
		  "  }\n");
    }
    if ($could_be_scalar) {	# Is a scalar allowed?
      $retstr .=  "  else\n"
	if $could_be_vec;
      $retstr .=  "    $argname = double(*$cvar);\n";
    }
    $retstr;
  };

$typemap_put_ptr{'complex < float >'} = sub { # Copy from intermediate array:
  my ($arg, $argname, $cvar) = @_; # Name the arguments.
  my $retstr = '';
  my $could_be_scalar = !@{$arg->{dimension}}; # Non-zero if this is definitely a vector.
  my $could_be_vec = @{$arg->{dimension}} || $arg->{vectorize}; # True if this might be a vector.
  if ($could_be_vec) {		# Could this be a vector?
    $retstr .=  "  if ($argname.IsArray())\n"
      if $could_be_scalar;	# Check (unless we know for sure).
    $retstr .= ("  {\n" .
		"    for (int _fidx = 0; _fidx < $argname.length(); ++_fidx, ++$cvar)\n" .
		"      $argname.ComplexPtr()[_fidx] = Tcomplex($cvar->real(), $cvar->imag());\n" .
		"  }\n");
  }
  if ($could_be_scalar) {	# Is a scalar allowed?
    $retstr .=  "  else\n"
      if $could_be_vec;
    $retstr .=  "    $argname = Tcomplex($cvar->real(), $cvar->imag());\n";
  }
  $retstr;
};

#
# Create an output argument vector of a given size.  Arguments:
# 1) The name of the argument.
# 2) The name of the C variable.
# 3) The name of a TDimPack structure.
# 4) The type.
#
$typemap_output_array_make{'double'} =
  sub { "    $_[0].rreserv($_[2]);\n    $_[1] = $_[0].RealPtr();\n" };

$typemap_output_array_make{'float'} = # Floats require an intermediate 
				# temporary vector.
  sub {
    my ($argname, $cvar, $dimpack, $type) = @_;	# Name the arguments.
    ("    $argname.rreserv($dimpack);\n" . # Save space for it.
     "    $cvar = (float *)alloca($dimpack.length() * sizeof (float));\n");
				# Make a temporary array of floats.
  };


$typemap_output_array_make{'int'} = $typemap_output_array_make{'unsigned'} =
  sub { "    $_[0].ireserv($_[2]);\n    $_[1] = ($_[3] *)$_[0].IntPtr();\n" };

$typemap_output_array_make{'complex < double >'} = sub {
  my ($argname, $cvar, $dimpack, $type) = @_;	# Name the arguments.
  ("    $argname.creserv($dimpack);\n" . # Save space for it.
   "    $cvar = ($type *)$argname.ComplexPtr();\n"); # Point to the space.
};

$typemap_output_array_make{'complex < float >'} = sub {
  my ($argname, $cvar, $dimpack, $type) = @_;	# Name the arguments.
  ("    $argname.creserv($dimpack);\n" . # Save space for it.
   "    $cvar = ($type *)alloca($dimpack.length() * sizeof ($type));\n");
				# Make a temporary array.
};


#
# Create an output scalar of a given size.  Arguments:
# 1) The name of the argument.
# 2) The name of the C variable.
# 3) The type.
#
$typemap_output_scalar_make{'double'} =
  sub { "    $_[0] = (double)0;\n    $_[1] = &$_[0].RealRef();\n" };
				# First use the operator= to indicate that it
				# is a scalar float, and then get a pointer
				# to its value.

$typemap_output_scalar_make{'float'} = # Floats require an intermediate.
  sub {
    my ($argname, $cvar, $type) = @_; # Name the arguments.
    ("    $cvar = (float *)alloca(sizeof (float));\n" .	# Make space for it.
     "    $argname.SetToVoid();\n"); # Flag that it's a scalar.
  };

$typemap_output_scalar_make{'int'} = $typemap_output_scalar_make{'unsigned'} =
  sub { "    $_[0] = 0;\n    $_[1] = ($_[2] *)&$_[0].IntRef();\n" };
				# Set the value to 0 so that we call Tela's
				# operator=(int) to set the type code.  Then
				# store the address.

$typemap_output_scalar_make{'char *'} = sub { # Nothing to do here.
};

$typemap_output_scalar_make{'char'} = sub { # Nothing to do here.
};

$typemap_output_scalar_make{'complex < double >'} = sub {
  my ($argname, $cvar, $type) = @_; # Name the arguments.
  ("    $argname = Tcomplex(0, 0);\n" .	# Flag it as complex.
   "    $cvar = ($type *)&$argname.ComplexRef();\n");
};

$typemap_output_scalar_make{'complex < float >'} = sub {
  my ($argname, $cvar, $type) = @_; # Name the arguments.
  ("    $argname = Tcomplex(0, 0);\n" .	# Flag it as complex.
   "    $cvar = ($type *)alloca(sizeof ($type));\n"); # Make temp. space.
};

#
# Below this point are the functions that get copied into the header of every
# module:
#
__DATA__
#ifndef __GNUC__
#include <alloca.h>
#endif
#include <complex.h>

//
// Function to get the length of a Tela object of arbitrary type.
// Arguments:
// 1) The tela object.
//
// If the object is a scalar, the array length is 1.
//
static inline int
_arraylen(const Tobject &obj) {
  if (obj.IsArray())
    return obj.length();
  else
    return 1;
}

//
// Get the n'th dimension of a Tela object.
// Arguments:
// 1) The Tela object.
// 2) Which dimension.  0 for the most quickly varying one, etc.
//    Note that because tela stores matrices in row major rather than
//    column major order, index 0 is actually the last index.
//
static int
_dim(const Tobject &obj, int idxno) {
  if (obj.HasStringFlag())	// Strings are always scalars.
    return 1;
  if (obj.IsArray() &&			// Do we actually have an array?
      obj.rank() > idxno)		// Has at least this # of dims?
    return obj.dims()[obj.rank() - idxno - 1]; // Return the size.
  else
    return 1;
}

//
// Return the number of dimensions in a Tela object.
// Arguments:
// 1) The object.
//
static int
_n_dims(const Tobject &obj)
{
  if (obj.HasStringFlag())
    return 0;			// Strings are always scalars.
  return (obj.IsArray() ? obj.rank() : 0);
}

//
// Load up an array (which may consist only of a single element) with
// the values from a Tela object.  Arguments:
// 1) The tela object.
// 2) Where to put the values.  The array should be allocated sufficiently
//    large so that all of the values in the tela object can be converted.
//
// These functions return 0 on error, or 1 if the operation succeded.
//
template <class FLOAT>
static int
_get_double(const Tobject &obj, FLOAT *arr) // float or double.
{
  int idx;
  switch (obj.kind())
  {
  case Kreal:			// Double precision:
    *arr = (FLOAT)obj.RealValue(); // Convert to correct type.
    break;
  case Kint:			// Integer:
    *arr = (FLOAT)obj.IntValue();
    break;

  case KRealArray:
    for (idx = 0; idx < _arraylen(obj); ++idx)
      *arr++ = (FLOAT)(obj.RealPtr()[idx]);
    break;
  case KIntArray:
    for (idx = 0; idx < _arraylen(obj); ++idx)
      *arr++ = (FLOAT)(obj.IntPtr()[idx]);
    break;

  default:
    return 0;			// Type conversion error.
  }

  return 1;			// No error.
}

//
// A version of the above function for complex numbers:
//
template <class FLOAT>
static int
_get_complex(const Tobject &obj, complex<FLOAT> *arr)
{
  int idx;
  switch (obj.kind())
  {
  case Kreal:			// Double precision:
    *arr = complex<FLOAT>(obj.RealValue());
    break;
  case Kint:			// Integer:
    *arr = complex<FLOAT>(obj.IntValue());
    break;

  case KRealArray:
    for (idx = 0; idx < _arraylen(obj); ++idx)
      *arr++ = complex<FLOAT>(obj.RealPtr()[idx]);
    break;
  case KIntArray:
    for (idx = 0; idx < _arraylen(obj); ++idx)
      *arr++ = complex<FLOAT>(obj.IntPtr()[idx]);
    break;

  case Kcomplex:		// Additions for complex numbers:
    *arr = complex<FLOAT>(obj.ComplexValue().real(),
			  obj.ComplexValue().imag());
    break;
  case KComplexArray:
    for (idx = 0; idx < _arraylen(obj); ++idx)
      *arr++ = complex<float>(obj.ComplexPtr()[idx].real(),
			      obj.ComplexPtr()[idx].imag());
    break;

  default:
    return 0;			// Type conversion error.
  }

  return 1;			// No error.
}

//
// Make sure that a Tela object is an integer, and return a pointer to it.
//
// Returns 0 if not an integer.
//
static int *
_get_int_ptr(const Tobject &obj)
{
  switch (obj.kind())
  {
  case Kint:
    return (int *)&((Tobject *)&obj)->IntRef();
  case KIntArray:
    return (int *)obj.IntPtr();
  default:
    return 0;			// Wrong type.
  }
}
