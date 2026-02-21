function write_orca_input(xyzfile,
                          charge,
                          multiplicity;
                          method="BP86 def2-TZVP",
                          solvent="CPCM(Water)",
                          maxcore=2000)

    open("input.inp","w") do io
        println(io, "%MaxCore $maxcore")
        println(io, "! $method $solvent")
        println(io, "* xyzfile $charge $multiplicity $xyzfile")
    end
end


function run_orca(xyzfile, charge, multiplicity;
                  method="BP86 def2-TZVP",
                  solvent="CPCM(Water)")

    write_orca_input(xyzfile, charge, multiplicity;
                     method=method,
                     solvent=solvent)

    run(`orca input.inp`)
end