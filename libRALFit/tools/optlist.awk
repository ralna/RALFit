# USAGE: ./optlist.awk options.90  | xclip
# options.f90 is a temporal files with the elements of the type nlls_options ...

awk '/^[[:space:]]*!/ {next}
     /^[[:space:]]*$/ {next}
     {print}' $1 | tr 'A-Z' 'a-z' | sed -e 's/kind[[:space:]]*=[[:space:]]*//' -e 's/[[:space:]]*=.*$//' \
     -e 's/^[[:space:]]*logical[[:space:]]*::[[:space:]]*\(.*\)/Write(adj,Fmt=99999) \"\1\"\nWrite(rec(),Fmt=99996) Adjustl(adj), options%\1/' \
     -e 's/^[[:space:]]*integer[[:space:]]*::[[:space:]]*\(.*\)/Write(adj,Fmt=99999) \"\1\"\nWrite(rec(),Fmt=99997) Adjustl(adj), options%\1/' \
     -e 's/^[[:space:]]*real[[:space:]]*(.*)[[:space:]]*::[[:space:]]*\(.*\)/Write(adj,Fmt=99999) \"\1\"\nWrite(rec(),Fmt=99998) Adjustl(adj), options%\1/'\
     | \
awk 'BEGIN {nrec=1}
     /rec/ {nrec=nrec+1; sub(/rec\(\)/,"rec("nrec")",$0); print "        " $0; next }
     {print "       ", $0}
     END {nrec = nrec + 1; print "        Write (rec("nrec"),Fmt=99994)"; print "        nrec =",nrec}'