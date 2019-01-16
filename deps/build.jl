using BinDeps

@BinDeps.setup
group = library_group("SatelliteScheduling")

libsofa_c = library_dependency("libsofa_c", aliases = ["libsofa_c", "libsofa"], runtime = true, group = group)

prefix = joinpath(BinDeps.depsdir(libsofa_c), "usr")
libdir = joinpath(BinDeps.depsdir(libsofa_c), "usr", "lib")
srcdir = joinpath(BinDeps.depsdir(libsofa_c), "src", "sofa")

provides(SimpleBuild,
    (@build_steps begin
    CreateDirectory(prefix)
    CreateDirectory(libdir)
        @build_steps begin
            ChangeDirectory(srcdir)
            `make clean`
            `make`
            `cp libsofa_c.so libsofa_c.dylib` 
            `mv libsofa_c.so $prefix/lib` 
            `mv libsofa_c.dylib $prefix/lib` 
            `make clean`
        end
    end), libsofa_c, os = :Unix)

@BinDeps.install Dict(:libsofa_c => :libsofa_c)