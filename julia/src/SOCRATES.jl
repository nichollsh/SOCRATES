"""
    SOCRATES

Julia interface to the SOCRATES radiative transfer codes

Calls C functions defined in `SOCRATES_C.f90` and `Str<name>_C.f90` Fortran modules
(using `ISO_C_BINDING` provided by Fortran 2003, with extensions in 2008).
"""
module SOCRATES

using OffsetArrays

const libSOCRATES_C = joinpath(@__DIR__, "../lib/libSOCRATES_C.so")

include("../gen/rad_pcf.jl")
include("../gen/gas_list_pcf.jl")
include("../gen/input_head_pcf.jl")

# Fortran TYPE wrappers
include("../gen/StrCtrl_JL.jl")
for member_name in [ # TYPE members of StrSpecData
        "Dim", "Basic", "Solar", "Rayleigh", "Gas", "Planck", "Cont", "ContGen", "Drop",
        "Aerosol", "Ice", "Var", "Photol", "Map",
    ]
    include("../gen/StrSpecData$(member_name)_JL.jl")
end
include("../gen/StrSpecData_JL.jl")
include("../gen/StrDim_JL.jl")
include("../gen/StrAtm_JL.jl")
include("../gen/StrCld_JL.jl")
include("../gen/StrAer_JL.jl")
include("../gen/StrBound_JL.jl")
include("../gen/StrOut_JL.jl")


##############################################
# StrCtrl
# radiance_core/def_control.F90
##############################################
 
"""
    StrCtrl
    
Opaque pointer to Fortran struct defined in
`radiance_core/def_control.F90`

Create and delete using:

    julia> control = SOCRATES.create_StrCtrl()
    julia> SOCRATES.delete_StrControl(control)

"""



function allocate_control(control::StrCtrl, sp::StrSpecData)
    ccall(
        (:PS_allocate_control, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid},), 
        control.cptr, sp.cptr,
    )
    return nothing
end

##############################################
# StrSpecData
# radiance_core/def_spectrum.F90
# interface_core/socrates_set_spectrum.F90
##############################################

"""
    StrSpecData

handle StrSpecData to the Fortran struct defined in 
`radiance_core/def_spectrum.F90`

"""

function read_spectrum(file_spectral::AbstractString, spectrum::StrSpecData)
    ccall(
        (:PS_read_spectrum, libSOCRATES_C),
        Cvoid, 
        (Cstring, Ptr{Cvoid},), 
        file_spectral, spectrum.cptr
    )
    return nothing
end

# interface_core/socrates_set_spectrum.F90
function set_spectrum(;
    n_instances::Union{Int, Ptr{Nothing}} = C_NULL,
    spectrum::Union{StrSpecData, Ptr{Nothing}} = C_NULL,
    spectrum_name::Union{AbstractString, Ptr{Nothing}} = C_NULL,
    spectral_file::Union{AbstractString, Ptr{Nothing}} = C_NULL,
    l_h2o      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_co2      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_o3       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_n2o      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_co       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch4      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_o2       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_no       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_so2      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_no2      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_nh3      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hno3     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_n2       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_cfc11    ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_cfc12    ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_cfc113   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hcfc22   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hfc125   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hfc134a  ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_cfc114   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_tio      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_vo       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_h2       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_he       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ocs      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_na       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_k        ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_feh      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_crh      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_li       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_rb       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_cs       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ph3      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_c2h2     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hcn      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_h2s      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ar       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_o        ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_n        ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_no3      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_n2o5     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hono     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ho2no2   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_h2o2     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_c2h6     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch3      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_h2co     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ho2      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hdo      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hcl      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hf       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_cosso    ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_tosso    ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_yosos    ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch3cho   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch3ooh   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch3coch3 ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch3cocho ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_chocho   ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_c2h5cho  ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_hoch2cho ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_c2h5coch3::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_mvk      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_macr     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_pan      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_ch3ono2  ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_sio      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_sio2     ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_fe       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_feo      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_na2      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_nao      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_mg       ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_mg2      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_mgo      ::Union{Bool, Ptr{Nothing}} = C_NULL,
    l_all_gasses::Union{Bool, Ptr{Nothing}} = C_NULL,
    wavelength_blue::Union{Float64, Ptr{Nothing}} = C_NULL,
)

    # NB: Fortran logical(c_bool) <-> C99 _Bool <->  Julia Cuchar

    ccall(
        (:PS_set_spectrum, libSOCRATES_C),
        Cvoid, 
        (
            Ref{Cint},
            Ptr{Cvoid},
            Cstring,
            Cstring,
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ref{Cuchar},
            Ptr{Cvoid}, # Ref{Cdouble} doesn't work (can't be set to C_NULL) ???
        ),
        n_instances,
        spectrum === C_NULL ? C_NULL : spectrum.cptr,
        spectrum_name,
        spectral_file,
        l_h2o        ,
        l_co2        ,
        l_o3         ,
        l_n2o        ,
        l_co         ,
        l_ch4        ,
        l_o2         ,
        l_no         ,
        l_so2        ,
        l_no2        ,
        l_nh3        ,
        l_hno3       ,
        l_n2         ,
        l_cfc11      ,
        l_cfc12      ,
        l_cfc113     ,
        l_hcfc22     ,
        l_hfc125     ,
        l_hfc134a    ,
        l_cfc114     ,
        l_tio        ,
        l_vo         ,
        l_h2         ,
        l_he         ,
        l_ocs        ,
        l_na         ,
        l_k          ,
        l_feh        ,
        l_crh        ,
        l_li         ,
        l_rb         ,
        l_cs         ,
        l_ph3        ,
        l_c2h2       ,
        l_hcn        ,
        l_h2s        ,
        l_ar         ,
        l_o          ,
        l_n          ,
        l_no3        ,
        l_n2o5       ,
        l_hono       ,
        l_ho2no2     ,
        l_h2o2       ,
        l_c2h6       ,
        l_ch3        ,
        l_h2co       ,
        l_ho2        ,
        l_hdo        ,
        l_hcl        ,
        l_hf         ,
        l_cosso      ,
        l_tosso      ,
        l_yosos      ,
        l_ch3cho     ,
        l_ch3ooh     ,
        l_ch3coch3   ,
        l_ch3cocho   ,
        l_chocho     ,
        l_c2h5cho    ,
        l_hoch2cho   ,
        l_c2h5coch3  ,
        l_mvk        ,
        l_macr       ,
        l_pan        ,
        l_ch3ono2    ,
        l_sio        ,
        l_sio2       ,
        l_fe         ,
        l_feo        ,
        l_na2        ,
        l_nao        ,
        l_mg         ,
        l_mg2        ,
        l_mgo        ,
        l_all_gasses,
        wavelength_blue === C_NULL ? C_NULL : Ref{Cdouble}(wavelength_blue),
    )

    return nothing
end

##############################################
# StrDim
# radiance_core/def_dimen.F90
##############################################
""" 
    StrDim

handle to the Fortran struct defined in 
`radiance_core/def_dimen.F90`.
"""

##############################################
# StrAtm
# radiance_core/def_atm.F90
##############################################


function allocate_atm(atm::StrAtm, dimen::StrDim, sp::StrSpecData)
    ccall(
        (:PS_allocate_atm, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},), 
        atm.cptr, dimen.cptr, sp.cptr,
    )
    return nothing
end

function deallocate_atm(atm::StrAtm)
    ccall(
        (:PS_deallocate_atm, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        atm.cptr,
    )
    return nothing
end

##############################################
# StrCld
# radiance_core/def_cld.F90
##############################################


function allocate_cld(cld::StrCld, dimen::StrDim, sp::StrSpecData)
    ccall(
        (:PS_allocate_cld, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},), 
        cld.cptr, dimen.cptr, sp.cptr,
    )
    return nothing
end

function deallocate_cld(cld::StrCld)
    ccall(
        (:PS_deallocate_cld, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        cld.cptr,
    )
    return nothing
end

function allocate_cld_prsc(cld::StrCld, dimen::StrDim, sp::StrSpecData)
    ccall(
        (:PS_allocate_cld_prsc, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},), 
        cld.cptr, dimen.cptr, sp.cptr,
    )
    return nothing
end

function deallocate_cld_prsc(cld::StrCld)
    ccall(
        (:PS_deallocate_cld_prsc, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        cld.cptr,
    )
    return nothing
end

##############################################
# StrAer
# radiance_core/def_aer.F90
##############################################


function allocate_aer(aer::StrAer, dimen::StrDim, sp::StrSpecData)
    ccall(
        (:PS_allocate_aer, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},), 
        aer.cptr, dimen.cptr, sp.cptr,
    )
    return nothing
end

function deallocate_aer(aer::StrAer)
    ccall(
        (:PS_deallocate_aer, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        aer.cptr,
    )
    return nothing
end

function allocate_aer_prsc(aer::StrAer, dimen::StrDim, sp::StrSpecData)
    ccall(
        (:PS_allocate_aer_prsc, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},), 
        aer.cptr, dimen.cptr, sp.cptr,
    )
    return nothing
end

function deallocate_aer_prsc(aer::StrAer)
    ccall(
        (:PS_deallocate_aer_prsc, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        aer.cptr,
    )
    return nothing
end

##############################################
# StrBound
# radiance_core/def_bound.F90
##############################################

function allocate_bound(bound::StrBound, dimen::StrDim, sp::StrSpecData)
    ccall(
        (:PS_allocate_bound, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},), 
        bound.cptr, dimen.cptr, sp.cptr,
    )
    return nothing
end

function deallocate_bound(bound::StrBound)
    ccall(
        (:PS_deallocate_bound, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        bound.cptr,
    )
    return nothing
end

##############################################
# StrOut
# radiance_core/def_out.F90
##############################################

function deallocate_out(radout::StrOut)
    ccall(
        (:PS_deallocate_out, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid},), 
        radout.cptr,
    )
    return nothing
end

###########################################
# calculate radiance
############################################

function radiance_calc(
    control::StrCtrl, dimen::StrDim, spectrum::StrSpecData, atm::StrAtm, cld::StrCld, aer::StrAer, bound::StrBound, radout::StrOut
)
    ccall(
        (:PS_radiance_calc, libSOCRATES_C),
        Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        control.cptr, dimen.cptr, spectrum.cptr, atm.cptr, cld.cptr, aer.cptr, bound.cptr, radout.cptr,
    )
    return nothing
end

###################################
# pretty printing
###################################

function dump_properties(io::IO, handle; sp="  ")
    for p in propertynames(handle)
        ps = getproperty(handle, p)
        if ps isa OffsetArray
            println(io, sp, "$p: $(size(ps))")
        else
            println(io, sp, "$p: $ps")
        end
    end
    return nothing
end

dump_properties(handle) = dump_properties(stdout, handle)

function dump_spectrum(spectrum)
    for p in propertynames(spectrum)
        ps = getproperty(spectrum, p)
        println("  $p: $ps")
        dump_properties(stdout, ps, sp="    ")        
    end    
    return nothing
end

######################################
# Test functions for argument passing
#####################################

function test_double_val(din::Float64)
    dout = ccall((:PS_test_double_val, libSOCRATES_C), Cdouble, (Cdouble, ), din)
    return dout
end

# din = C_NULL -> Fortran optional argument not present
function test_double_ref(din::Union{Float64, Ptr{Nothing}} = C_NULL)
    # dout = ccall((:PS_test_double_ref, libSOCRATES_C), Cdouble, (Ref{Cdouble}, ), din)

    # TODO -  Ref{Cdouble} apparently can't be set to C_NULL, so 
    # we have to pass the argument as Ptr{Cvoid}
    # cf Ref{Cuchar} (Fortran c_bool <-> C99 _Bool ) which can be set to C_NULL
    din_C = din === C_NULL ? C_NULL : Ref{Cdouble}(din)
    dout = ccall((:PS_test_double_ref, libSOCRATES_C), Cdouble, (Ptr{Cvoid}, ), din_C)
    return dout
end

end
