while (<>) {
  s/(^|[ \t])([^-][^ ]*[.]a)/ -Wl,$2/g;
  print;
}
