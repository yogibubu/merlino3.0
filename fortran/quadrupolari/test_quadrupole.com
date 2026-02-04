***, NH3FH -- EFG tensor and quadrupole coupling for 14N
memory,4000,m
wf,charge=0,spin=0

geometry={
6

 N    -1.371091    0.000000   -0.000000
 H    -1.744143   -0.471455    0.816584
 H    -1.744143   -0.471455   -0.816584
 H    -1.744143    0.942911   -0.000000
 H     0.315211    0.000000    0.000000
 F     1.271432    0.000000    0.000000
}

basis=cc-pVTZ

{hf; expec,fg}
{ccsd(t); expec,fg}

