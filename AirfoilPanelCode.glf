package require PWI_Glyph

proc determ {arrayVarName n} {

    upvar $arrayVarName arr

    for {set i 1} {$i <= $n} {incr i} {

        for {set j 1} {$j <= $n} {incr j} {

            set a($i,$j) $arr($i,$j)

        }

    }

    set m 1
    set k [expr {$m+1}]
   
    while {$m <= $n} {

        if {$m == [expr {$n-1}]} {

            set determ 1.0

        }
 
        for {set i $k} {$i <= $n} {incr i} {

            set ratio [expr {$a($i,$m)/$a($m,$m)}]

            for {set j $k} {$j <= $n} {incr j} {

                set a($i,$j) [expr {$a($i,$j)-$a($m,$j)*$ratio}]

            }

        }

        incr m
        set k [expr {$m+1}]

    }

    for {set l 1} {$l <= $n} {incr l} {

        set determ [expr {$determ*$a($l,$l)}]

    }

    return $determ

}


proc cramer {cVarName aVarName xVarName n} {

    upvar $cVarName c
    upvar $aVarName a
    upvar $xVarName x

    set denom [determ c $n]    

    for {set k 1} {$k <= $n} {incr k} {

        for {set i 1} {$i <= $n} {incr i} {

            for {set j 1} {$j <= $n} {incr j} {

                set cc($i,$j) $c($i,$j)

            }

        }

        for {set i 1} {$i <= $n} {incr i} {
        
            set cc($i,$k) $a($i)

        }

        set x($k) [expr {[determ cc $n]/$denom}]

    }

}


proc naca4gen {naca n xVarName yVarName} {

    upvar $xVarName xb
    upvar $yVarName yb

	# AIRFOIL INPUTS
	# -----------------------------------------------
	set m [expr {[string index $naca 0]/100.0}]; # maximum camber 
	set p [expr {[string index $naca 1]/10.0}] ; # maximum camber location
	set a [string index $naca 2]
	set b [string index $naca 3]
	set c "$a$b"
	set t [expr {$c/100.0}]                    ; # maximum thickness
	
	# GENERATE AIRFOIL COORDINATES
	# -----------------------------------------------
	# Airfoil step size
	set ds [expr {1.0/$n}]
	
	# Check if airfoil is symmetric or cambered
	if {$m == 0 && $p == 0 || $m == 0 || $p == 0} {set symm 1} else {set symm 0}
	
	# Get x coordinates
	for {set i 0} {$i < [expr {1+$ds}]} {set i [expr {$i+$ds}]} {lappend x $i}
	
	# Calculate mean camber line and thickness distribution
    set yc 0
	foreach xx $x {
	
		# Mean camber line definition for symmetric geometry
		if {$symm == 1} {lappend yc 0}
	
		# Mean camber line definition for cambered geometry
		if {$symm == 0 && $xx <= $p} {
			lappend yc [expr {($m/($p**2))*(2*$p*$xx-$xx**2)}]
		} elseif {$symm == 0 && $xx > $p} {
			lappend yc [expr {($m/((1-$p)**2)*(1-2*$p+2*$p*$xx-$xx**2))}]
		}
	
		# Thickness distribution
		lappend yt [expr {($t/0.20)*(0.29690*sqrt($xx)-0.12600*$xx- \
		                  0.35160*$xx**2+0.28430*$xx**3-0.1036*$xx**4)}]
	
		# Theta
		set dy [expr {[lindex $yc end] - [lindex $yc end-1]}]
		set th [expr {atan($dy/$ds)}]
	
		# Upper x and y coordinates
		lappend xu [expr {$xx-[lindex $yt end]*sin($th)}]
		lappend yu [expr {[lindex $yc end]+[lindex $yt end]*cos($th)}]
	
		# Lower x and y coordinates
		lappend xl [expr {$xx+[lindex $yt end]*sin($th)}]
		lappend yl [expr {[lindex $yc end]-[lindex $yt end]*cos($th)}]
	
	}

    # Create list of all x and y airfoil coordinates
    set xs [concat [lreverse $xl] [lrange $xu 1 end]]
    set ys [concat [lreverse $yl] [lrange $yu 1 end]]
    
    # Construct panel boundary point arrays
    set i 1
    foreach x $xs y $ys {

        set xb($i) $x
        set yb($i) $y
        incr i 

    }

}


set airfoil    2412; # NACA 4-series airfoil
set halfPanels 20  ; # half the number of panels
set alphaDeg   8.0 ; # angle of attack

# coordinates start at the trailing edge and march counter clockwise
# and the trailing edge is included twice

naca4gen $airfoil $halfPanels xb yb

set m     [expr {[array size xb]-1}]
set mp1   [expr {$m+1}]
set pi    [expr {4.0*atan(1.0)}]
set alpha [expr {$alphaDeg*$pi/180.0}]

for {set i 1} {$i <= $m} {incr i} {

    set ip1 [expr {$i+1}]

    set x($i) [expr {0.5*($xb($i)+$xb($ip1))}]
    set y($i) [expr {0.5*($yb($i)+$yb($ip1))}]

    set dx($i) [expr {$xb($ip1)-$xb($i)}]
    set dy($i) [expr {$yb($ip1)-$yb($i)}]
    
    set s($i)      [expr {sqrt(pow($xb($ip1)-$xb($i),2)+pow($yb($ip1)-$yb($i),2))}]
    set theta($i)  [expr {atan2($yb($ip1)-$yb($i),$xb($ip1)-$xb($i))}]
    set sine($i)   [expr {sin($theta($i))}]
    set cosine($i) [expr {cos($theta($i))}]
    set rhs($i)    [expr {sin($theta($i)-$alpha)}]

}

for {set i 1} {$i <= $m} {incr i} {

    for {set j 1} {$j <= $m} {incr j} {

        if {$i == $j} {

            set cn1($i,$j) -1.0
            set cn2($i,$j) 1.0
            set ct1($i,$j) [expr {0.5*$pi}]
            set ct2($i,$j) [expr {0.5*$pi}]

        } else {

            set a [expr {-1.0*($x($i)-$xb($j))*$cosine($j)-($y($i)-$yb($j))*$sine($j)}]
            set b [expr {pow($x($i)-$xb($j),2)+pow($y($i)-$yb($j),2)}]
            set c [expr {sin($theta($i)-$theta($j))}]
            set d [expr {cos($theta($i)-$theta($j))}]
            set e [expr {($x($i)-$xb($j))*$sine($j)-($y($i)-$yb($j))*$cosine($j)}]
            set f [expr {log(1.0+$s($j)*($s($j)+2.0*$a)/$b)}]
            set g [expr {atan2($e*$s($j),$b+$a*$s($j))}]
            set p [expr {($x($i)-$xb($j))*sin($theta($i)-2.0*$theta($j))\
                        +($y($i)-$yb($j))*cos($theta($i)-2.0*$theta($j))}]
            set q [expr {($x($i)-$xb($j))*cos($theta($i)-2.0*$theta($j))\
                        -($y($i)-$yb($j))*sin($theta($i)-2.0*$theta($j))}]

            set cn2($i,$j) [expr {$d+0.5*$q*$f/$s($j)-($a*$c+$d*$e)*$g/$s($j)}]
            set cn1($i,$j) [expr {0.5*$d*$f+$c*$g-$cn2($i,$j)}]

            set ct2($i,$j) [expr {$c+0.5*$p*$f/$s($j)+($a*$d-$c*$e)*$g/$s($j)}]
            set ct1($i,$j) [expr {0.5*$c*$f-$d*$g-$ct2($i,$j)}] 

        }

    }

}

for {set i 1} {$i <= $m} {incr i} {

    set an($i,1)    $cn1($i,1)
    set an($i,$mp1) $cn2($i,$m)

    set at($i,1)    $ct1($i,1)
    set at($i,$mp1) $ct2($i,$m)

    for {set j 2} {$j <= $m} {incr j} {
    
        set an($i,$j) [expr {$cn1($i,$j)+$cn2($i,[expr {$j-1}])}]
        set at($i,$j) [expr {$ct1($i,$j)+$ct2($i,[expr {$j-1}])}]

    }

}

set an($mp1,1)    1.0
set an($mp1,$mp1) 1.0

for {set j 2} {$j <= $m} {incr j} {

    set an($mp1,$j) 0.0

}

set rhs($mp1) 0.0

array set gama {}
cramer an rhs gama $mp1

set cfx 0
set cfy 0
for {set i 1} {$i <= $m} {incr i} {

    set v($i) [expr {cos($theta($i)-$alpha)}]

    for {set j 1} {$j <= $mp1} {incr j} {

        set v($i)  [expr {$v($i)+$at($i,$j)*$gama($j)}]
        set cp($i) [expr {1.0-pow($v($i),2)}]

    }

    puts $cp($i)

    set cfx [expr {$cfx+$cp($i)*$dy($i)}]
    set cfy [expr {$cfy-$cp($i)*$dx($i)}]

}

set CL [expr {$cfy*cos($alpha)-$cfx*sin($alpha)}]
puts ""
puts "Coefficient of Lift: $CL"
