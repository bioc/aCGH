.First.lib <-
    function(libname, pkgname, where)
    library.dynam("aCGH", pkgname, libname)

.Last.lib <-
    function(libpath)
    dyn.unload(file.path(libpath,
                         "libs",
                         paste("aCGH",
                               .Platform$"dynlib.ext",
                               sep = "")))
