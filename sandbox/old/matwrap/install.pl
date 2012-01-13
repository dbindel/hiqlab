#!/usr/local/bin/perl
#
# This script installs the bib_gen program in all the appropriate places
# after changing the paths in the scripts.  You probably shouldn't run this
# script directly; see the instructions in the file INSTALL for how to use it.
#
# Usage:
#	install.pl [-prefix dir] list-of-files
#
# -prefix specifies the output directory, e.g., /usr/local.  If this is
# blank, then we install the product in place.
#
# This script is intended to be called from a makefile, so don't use it
# directly.
#
require 5.002;
use Config;
use Cwd;
use FileHandle;

$prefix = '/usr/local';		# Set up default values of variables.

#
# First find out where perl is in our path:
#
$perlbin = "$Config{'bin'}/perl";
unless (-x $perlbin) {		# Not correctly specified by Config?
  foreach (split(/:/, $ENV{'PATH'})) { # Look in the path.
    if (-x "$_/perl") {		# Found it?
      $perlbin = "$_/perl";
      last;
    }
  }
}

-x $perlbin or die("Can't find executable perl from Config or from path\n");

#
# Now find out where this is being installed from:
#
($srcdir = $0) =~ s([^/]+$)();	# Get the directory name.
unless ($srcdir =~ m(^/)) {	# Is the full path specified?
  $srcdir = cwd . ($srcdir =~ /^\.?/ ? "/$srcdir" : ""); # Make the path absolute.
  $srcdir =~ s@/\.?$@@;		# Trim useless /. off the end.
}

#
# Parse command line options:
#
@files = ();
while ($_ = shift(@ARGV)) {	# Look at the next argument.
  if (/^-/) {			# Is this an option?
    $_ eq '-prefix' or die("usage: install.pl [-prefix dir] list-of-files\n");
    $prefix = shift(@ARGV);	# Take the next argument as the prefix.
    $prefix =~ s(/$)();		# Strip off a trailing slash, if present.
  } else {			# It's a file:
    push(@files, $_);		# Just remember it.
  }
}

unless ($prefix =~ m(^/)) {	# Is the full path specified?
  $prefix = cwd . ($prefix =~ /^\.?$/ ? "/$prefix" : ""); # Make the path absolute.
  $prefix =~ s@/\.$@@;		# Trim useless /. off the end.
}

if ($prefix) {			# Was a new location specified?
  -d $prefix || mkdir($prefix, 0755) ||
    die("Could not create $prefix--$!\n");
  -d "$prefix/man" || mkdir("$prefix/man", 0755) ||
    die("Could not create $prefix/man--$!\n");
  -d "$prefix/man/man1" || mkdir("$prefix/man/man1", 0755) ||
    die("Could not create $prefix/man--$!\n");
} else {
  $prefix = $srcdir;		# Prefix not specified--leave files in place.
}

#
# Start installing the files:
#
foreach (@files) {
  &copy_with_substitutions($_);
}

#
# Now generate the man pages:
#
foreach (grep /\.pod$/, @files) { # Take all of the pod files:
  my $manfile = $_;
  $manfile =~ s(.*/)();		# Strip off the directory name.
  $manfile =~ s/\.pod$/.1/;	# Replace the suffix.
  system("pod2man --section=1 $srcdir/$_ > $prefix/man/man1/$manfile");
}

#
# Copy a file, making appropriate substitutions.  Since none of the files
# we copy is large, we just read it into memory and then write it out.
# Arguments:
# 1) The name.  It is copied from $srcdir to $prefix.
#
sub copy_with_substitutions {
  my ($fname) = @_;		# Access the arguments.

  my $mode = (stat("$srcdir/$fname"))[2]; # Get the file protections.

  if ($srcdir eq $prefix) {	# Overwriting the output?
    rename("$srcdir/$fname", "$srcdir/$fname~"); # For safety, keep the old one around.
    open(INFILE, "$srcdir/$fname~") ||
      die("Can't open file $srcdir/$fname~--$!\n");
  } else {			# Copying to a different place:
#
# Make sure the directory exists:
#
    my @dirs = split(/\//, $fname); # Get the directories.
    pop(@dirs);			# Get rid of the file name at the end.
    my $dirstr = $prefix;	# Start with the directories that already exist.
    foreach (@dirs) {
      -d "$dirstr/$_" || mkdir("$dirstr/$_", 0755) || # Make the directory.
	die("Could not create directory $dirstr/$_--$!\n");
      $dirstr .= "/$_";		# Accumulate the directory string.
    }

    open(INFILE, "$srcdir/$fname") ||
      die("Can't open file $srcdir/$fname--$!\n");
  }

  local $/ = undef;		# Slurp in the whole file.
  $_ = <INFILE>;
  close(INFILE);		# Done with this file.

  s/^#!.*perl/#!$perlbin/;	# Change the path to perl.
  s(\$PREFIX = \"[^\"]+\")(\$PREFIX = \"$prefix\");
				# Change the path to the installed files.

  open(OUTFILE, "> $prefix/$fname") ||
    die("Can't create $prefix/$fname--$!\n");
  print OUTFILE $_;		# Print the modified text.
  close(OUTFILE);		# Done with this file.

  chmod $mode, "$prefix/$fname"; # Make sure it has the correct protections.
}
