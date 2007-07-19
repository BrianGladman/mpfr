# convert from Gonnet's FPAccuracy data sets to mpfr format
# http://www.inf.ethz.ch/personal/gonnet/FPAccuracy/all.tar.Z

# 1 - cut the lines from (say) C/acos.c, and remove the eps field
#     (hint: cut -d" " -f1,2,4,5 /tmp/acos.c > /tmp/acos2.c)
# 2 - edit the infile and outfile lines below, and run
#     maple -q < gonnet.mpl 

infile := "/tmp/acos2.c":
outfile := "acos":

###################### don't edit below this line #############################

foo := proc(arg_m, val_m, arg_e, val_e, fp)
   fprintf (fp, "53 53 n ", 53);
   to_hex(arg_m, arg_e, fp);
   fprintf (fp, " ");
   to_hex(val_m, -val_e, fp);
   fprintf (fp, "\n");
end:

to_hex := proc(m, e, fp)
   if m<0 then fprintf (fp, "-") fi;
   fprintf (fp, "0x%sp%d", convert(abs(m),hex), e);
end:

fp := fopen (outfile, WRITE):

l := readdata(infile, integer, 4):
for e in l do foo(op(e), outfile) od:

fclose (fp);

quit;

