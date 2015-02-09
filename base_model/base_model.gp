set term x11 0

set multiplot layout 1,2 title "ok"

set title "s_i"
plot 'state' using 0:1 with lines

set title "a"
plot 'state' using 0:2 with lines

unset multiplot
 
pause -1

set multiplot layout 2,2 title "ok"

set title "r"
plot 'state' using 0:3 with lines

set title "e_t"
plot 'state' using 0:4 with lines

set title "e_m"
plot 'state' using 0:5 with lines

set title "q"
plot 'state' using 0:6 with lines

unset multiplot

pause -1

set multiplot layout 2,2 title "ok"

set title "m_r"
plot 'state' using 0:7 with lines

set title "m_t"
plot 'state' using 0:8 with lines

set title "m_m"
plot 'state' using 0:9 with lines

set title "m_q"
plot 'state' using 0:10 with lines

unset multiplot

pause -1

set multiplot layout 2,2 title "ok"

set title "c_r"
plot 'state' using 0:11 with lines

set title "c_t"
plot 'state' using 0:12 with lines

set title "c_m"
plot 'state' using 0:13 with lines

set title "c_q"
plot 'state' using 0:14 with lines

unset multiplot

pause -1