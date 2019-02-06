#!/bin/sh

echo $2 | split -b60
optstr=`cat xaa`
rm xa*

cmp=$1
opt1=$optstr
opt2=""
opt3=""
host=`hostname`
commit=`git rev-parse HEAD `
set `date`
cat > date.h <<EOF
      character(80), parameter :: vers =  &
     &  "Version compiled $3 $2 $6 at $4 on "// &
     &  "$host "
      character(100), parameter :: cmp =  &
     &  "$cmp "
      character(60), parameter :: opt1 = &
     &  "$opt1 "
      character(60), parameter :: opt2 = &
     &  "$opt2 "
      character(60), parameter :: opt3 = &
     &  "$opt3 "
      character(40), parameter :: commit = &
     &  "$commit"

EOF


