#!/bin/tcsh 

set copystmt = $1
shift

foreach filename ($*)
  cp $filename $filename-bak
  cat $copystmt $filename-bak > $filename
end
