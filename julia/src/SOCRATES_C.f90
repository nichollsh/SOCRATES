! C interface to SOCRATES radiative transfer codes.
!
! Uses ISO_C_BINDING provided by Fortran 2003, 2008.
!
! Overall approach:
!   1) Fortran structs (StrCtrl, StrSpecData, etc) are passed as opaque C void* pointers,
!     which are created and deleted, and accessed / modified using C functions (create_StrCtrl,
!     delete_StrCtrl, StrCtrl_get_... , StrCtrl_set_...) provided by auto-generated wrappers.
!
!   2) The SOCRATES 'interface to the calling model' API as far as possible is a thin C wrapper around 
!     the corresponding Fortran functions in SOCRATES src/radiance_core. These are hand-generated.
!
module SOCRATES_C

use, intrinsic :: iso_c_binding

USE UTILITIES_CF

USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
USE def_spectrum
USE socrates_set_spectrum, only: set_spectrum
USE def_dimen
USE def_atm
USE def_cld
USE def_aer
USE def_bound
USE def_out

implicit none

private



contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrCtrl      
    !   control options initially read in via a namelist
    ! radiance_core/def_control.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! C-bound wrapper for allocate_control
    subroutine PS_allocate_control(control_Cptr, spectrum_Cptr) bind(C, name='PS_allocate_control')
        implicit none

        type(c_ptr), value, intent(in)          :: control_Cptr
        type(c_ptr), value, intent(in)          :: spectrum_Cptr

        ! local variables
        type(StrCtrl), pointer                  :: control
        type(StrSpecData), pointer              :: spectrum
        
        ! convert types
        call C_F_POINTER(control_Cptr, control)
        call C_F_POINTER(spectrum_Cptr, spectrum)
                
        call allocate_control(control, spectrum)

    end subroutine PS_allocate_control
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrSpecData
    !   spectral discretisation and optical properties read in from the spectral file
    ! radiance_core/def_spectrum.F90
    ! interface_core/socrates_set_spectrum.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    ! C-bound wrapper for radiance_core/read_spectrum
    subroutine PS_read_spectrum(file_spectral_Cstr, Sp_Cptr) bind(C, name='PS_read_spectrum')
        implicit none

        type(c_ptr), value, intent(in)          :: file_spectral_Cstr
        type(c_ptr), value                      :: Sp_Cptr

        ! local variables
        type(StrSpecData), pointer              :: Sp
        character(len=:), allocatable           :: file_spectral

        ! convert types
        call C_F_POINTER(Sp_Cptr, Sp)
        file_spectral = c_to_f_string(file_spectral_Cstr)
        
        call read_spectrum(file_spectral, Sp)

    end subroutine PS_read_spectrum

    ! C-bound wrapper for interface_core/socrates_set_spectrum set_spectrum
    subroutine PS_set_spectrum(n_instances, spectrum_Cptr, spectrum_name_Cstr, spectral_file_Cstr, &
        l_h2o, l_co2, l_o3, l_n2o, l_co, l_ch4, l_o2, l_no, l_so2, l_no2, l_nh3, &
        l_hno3, l_n2, l_cfc11, l_cfc12, l_cfc113, l_hcfc22, l_hfc125, l_hfc134a, &
        l_cfc114, l_tio, l_vo, l_h2, l_he, l_ocs, l_na, l_k, l_feh, l_crh, l_li, &
        l_rb, l_cs, l_ph3, l_c2h2, l_hcn, l_h2s, l_ar, l_o, l_n, l_no3, l_n2o5, &
        l_hono, l_ho2no2, l_h2o2, l_c2h6, l_ch3, l_h2co, l_ho2, l_hdo, l_hcl, &
        l_hf, l_cosso, l_tosso, l_yosos, l_ch3cho, l_ch3ooh, l_ch3coch3, &
        l_ch3cocho, l_chocho, l_c2h5cho, l_hoch2cho, l_c2h5coch3, l_mvk, l_macr, &
        l_pan, l_ch3ono2, l_sio, l_sio2, l_fe, l_feo, l_na2, l_nao, l_mg, l_mg2, l_mgo, l_all_gases, &
        wavelength_blue) bind(C, name='PS_set_spectrum')   
    
        implicit none
            
        ! Number of instances of the spectrum type (to allocate spectrum_array)
        integer(c_int), intent(in), optional :: n_instances
        
        ! Spectral data:
        type(c_ptr), value, intent(in)              :: spectrum_Cptr
       
        type(c_ptr), value, intent(in)              :: spectrum_name_Cstr
        type(c_ptr), value, intent(in)              :: spectral_file_Cstr
        
        logical(c_bool), intent(in), optional :: &
            l_h2o, l_co2, l_o3, l_n2o, l_co, l_ch4, l_o2, l_no, l_so2, l_no2, l_nh3, &
            l_hno3, l_n2, l_cfc11, l_cfc12, l_cfc113, l_hcfc22, l_hfc125, l_hfc134a, &
            l_cfc114, l_tio, l_vo, l_h2, l_he, l_ocs, l_na, l_k, l_feh, l_crh, l_li, &
            l_rb, l_cs, l_ph3, l_c2h2, l_hcn, l_h2s, l_ar, l_o, l_n, l_no3, l_n2o5, &
            l_hono, l_ho2no2, l_h2o2, l_c2h6, l_ch3, l_h2co, l_ho2, l_hdo, l_hcl, &
            l_hf, l_cosso, l_tosso, l_yosos, l_ch3cho, l_ch3ooh, l_ch3coch3, &
            l_ch3cocho, l_chocho, l_c2h5cho, l_hoch2cho, l_c2h5coch3, l_mvk, l_macr, &
            l_pan, l_ch3ono2, l_sio, l_sio2, l_fe, l_feo, l_na2, l_nao, l_mg, l_mg2, l_mgo, & 
            l_all_gases
        
        real(c_double), intent(in), optional :: wavelength_blue

        ! local variables
        type(StrSpecData), pointer              :: spectrum => null()
        character(len=:), allocatable           :: spectrum_name
        character(len=:), allocatable           :: spectral_file

        ! convert types
        ! TODO check: null pointer is equivalent to optional argument not present in Fortran 2008 ?
        if (c_associated(spectrum_Cptr)) then
            call C_F_POINTER(spectrum_Cptr, spectrum)
        end if
        ! variable not allocated is equivalent to optional argument not present in Fortran 2008
        if (c_associated(spectrum_name_Cstr)) then 
            spectrum_name = c_to_f_string(spectrum_name_Cstr)
        end if
        if (c_associated(spectral_file_Cstr)) then
            spectral_file = c_to_f_string(spectral_file_Cstr)
        end if

        ! Read spectral file
        CALL set_spectrum( &
            n_instances = n_instances, &
            spectrum      = spectrum, &
            spectral_file = spectral_file, &
            spectrum_name = spectrum_name, &
            l_h2o         = logical(l_h2o      ), & ! c_bool is not the same as default 'kind' of Fortran logical 
            l_co2         = logical(l_co2      ), & 
            l_o3          = logical(l_o3       ), & 
            l_n2o         = logical(l_n2o      ), & 
            l_co          = logical(l_co       ), & 
            l_ch4         = logical(l_ch4      ), & 
            l_o2          = logical(l_o2       ), & 
            l_no          = logical(l_no       ), & 
            l_so2         = logical(l_so2      ), & 
            l_no2         = logical(l_no2      ), & 
            l_nh3         = logical(l_nh3      ), & 
            l_hno3        = logical(l_hno3     ), & 
            l_n2          = logical(l_n2       ), & 
            l_cfc11       = logical(l_cfc11    ), & 
            l_cfc12       = logical(l_cfc12    ), & 
            l_cfc113      = logical(l_cfc113   ), & 
            l_hcfc22      = logical(l_hcfc22   ), & 
            l_hfc125      = logical(l_hfc125   ), & 
            l_hfc134a     = logical(l_hfc134a  ), & 
            l_cfc114      = logical(l_cfc114   ), & 
            l_tio         = logical(l_tio      ), & 
            l_vo          = logical(l_vo       ), & 
            l_h2          = logical(l_h2       ), & 
            l_he          = logical(l_he       ), & 
            l_ocs         = logical(l_ocs      ), & 
            l_na          = logical(l_na       ), & 
            l_k           = logical(l_k        ), & 
            l_feh         = logical(l_feh      ), & 
            l_crh         = logical(l_crh      ), & 
            l_li          = logical(l_li       ), & 
            l_rb          = logical(l_rb       ), & 
            l_cs          = logical(l_cs       ), & 
            l_ph3         = logical(l_ph3      ), & 
            l_c2h2        = logical(l_c2h2     ), & 
            l_hcn         = logical(l_hcn      ), & 
            l_h2s         = logical(l_h2s      ), & 
            l_ar          = logical(l_ar       ), & 
            l_o           = logical(l_o        ), & 
            l_n           = logical(l_n        ), & 
            l_no3         = logical(l_no3      ), & 
            l_n2o5        = logical(l_n2o5     ), & 
            l_hono        = logical(l_hono     ), & 
            l_ho2no2      = logical(l_ho2no2   ), & 
            l_h2o2        = logical(l_h2o2     ), & 
            l_c2h6        = logical(l_c2h6     ), & 
            l_ch3         = logical(l_ch3      ), & 
            l_h2co        = logical(l_h2co     ), & 
            l_ho2         = logical(l_ho2      ), & 
            l_hdo         = logical(l_hdo      ), & 
            l_hcl         = logical(l_hcl      ), & 
            l_hf          = logical(l_hf       ), & 
            l_cosso       = logical(l_cosso    ), &
            l_tosso       = logical(l_tosso    ), &
            l_yosos       = logical(l_yosos    ), &
            l_ch3cho      = logical(l_ch3cho   ), &
            l_ch3ooh      = logical(l_ch3ooh   ), &
            l_ch3coch3    = logical(l_ch3coch3 ), &
            l_ch3cocho    = logical(l_ch3cocho ), &
            l_chocho      = logical(l_chocho   ), &
            l_c2h5cho     = logical(l_c2h5cho  ), &
            l_hoch2cho    = logical(l_hoch2cho ), &
            l_c2h5coch3   = logical(l_c2h5coch3), &
            l_mvk         = logical(l_mvk      ), &
            l_macr        = logical(l_macr     ), &
            l_pan         = logical(l_pan      ), &
            l_ch3ono2     = logical(l_ch3ono2  ), &
            l_sio         = logical(l_sio      ), &
            l_sio2        = logical(l_sio2     ), &
            l_fe          = logical(l_fe       ), & 
            l_feo         = logical(l_feo      ), & 
            l_na2         = logical(l_na2      ), & 
            l_nao         = logical(l_nao      ), & 
            l_mg          = logical(l_mg       ), & 
            l_mg2         = logical(l_mg2      ), & 
            l_mgo         = logical(l_mgo      ), & 
            l_all_gases   = logical(l_all_gases), &
            wavelength_blue = wavelength_blue )

    end subroutine PS_set_spectrum


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrDim
    !   dimensions for arrays
    ! radiance_core/def_dimen.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrAtm
    !   grid discretisation and atmospheric profiles of thermodynamic quantities and gas amounts
    ! radiance_core/def_atm.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! C-bound wrapper for allocate_atm
    subroutine PS_allocate_atm(atm_Cptr, dimen_Cptr, sp_Cptr) bind(C, name='PS_allocate_atm')
        implicit none

        type(c_ptr), value, intent(in)          :: atm_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: sp_Cptr

        ! local variables
        type(StrAtm), pointer                   :: atm
        type(StrDim), pointer                 :: dimen
        type(StrSpecData), pointer              :: sp
        
        ! convert types
        call C_F_POINTER(atm_Cptr, atm)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(sp_Cptr, sp)
                
        call allocate_atm(atm, dimen, sp)

    end subroutine PS_allocate_atm

    ! C-bound wrapper for deallocate_atm
    subroutine PS_deallocate_atm(atm_Cptr) bind(C, name='PS_deallocate_atm')
        implicit none

        type(c_ptr), value, intent(in)          :: atm_Cptr

        ! local variables
        type(StrAtm), pointer                   :: atm
        
        ! convert types
        call C_F_POINTER(atm_Cptr, atm)
                
        call deallocate_atm(atm)

    end subroutine PS_deallocate_atm
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrCld
    !   cloud fields (fractions, mixing ratios and sub-grid structure)
    ! radiance_core/def_cld.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! C-bound wrapper for allocate_cld
    subroutine PS_allocate_cld(cld_Cptr, dimen_Cptr, sp_Cptr) bind(C, name='PS_allocate_cld')
        implicit none

        type(c_ptr), value, intent(in)          :: cld_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: sp_Cptr

        ! local variables
        type(StrCld), pointer                   :: cld
        type(StrDim), pointer                 :: dimen
        type(StrSpecData), pointer              :: sp
        
        ! convert types
        call C_F_POINTER(cld_Cptr, cld)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(sp_Cptr, sp)
                
        call allocate_cld(cld, dimen, sp)

    end subroutine PS_allocate_cld

    ! C-bound wrapper for deallocate_cld
    subroutine PS_deallocate_cld(cld_Cptr) bind(C, name='PS_deallocate_cld')
        implicit none

        type(c_ptr), value, intent(in)          :: cld_Cptr

        ! local variables
        type(StrCld), pointer                   :: cld
        
        ! convert types
        call C_F_POINTER(cld_Cptr, cld)
                
        call deallocate_cld(cld)

    end subroutine PS_deallocate_cld
    
    ! C-bound wrapper for allocate_cld_prsc
    subroutine PS_allocate_cld_prsc(cld_Cptr, dimen_Cptr, sp_Cptr) bind(C, name='PS_allocate_cld_prsc')
        implicit none

        type(c_ptr), value, intent(in)          :: cld_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: sp_Cptr

        ! local variables
        type(StrCld), pointer                   :: cld
        type(StrDim), pointer                 :: dimen
        type(StrSpecData), pointer              :: sp
        
        ! convert types
        call C_F_POINTER(cld_Cptr, cld)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(sp_Cptr, sp)
                
        call allocate_cld_prsc(cld, dimen, sp)

    end subroutine PS_allocate_cld_prsc

    ! C-bound wrapper for deallocate_cld_prsc
    subroutine PS_deallocate_cld_prsc(cld_Cptr) bind(C, name='PS_deallocate_cld_prsc')
        implicit none

        type(c_ptr), value, intent(in)          :: cld_Cptr

        ! local variables
        type(StrCld), pointer                   :: cld
        
        ! convert types
        call C_F_POINTER(cld_Cptr, cld)
                
        call deallocate_cld_prsc(cld)

    end subroutine PS_deallocate_cld_prsc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrAer
    !   aerosol fields (mixing ratios for CLASSIC aerosols,
    !   optical properties for GLOMAPMODE aerosols)
    ! radiance_core/def_aer.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! C-bound wrapper for allocate_aer
    subroutine PS_allocate_aer(aer_Cptr, dimen_Cptr, sp_Cptr) bind(C, name='PS_allocate_aer')
        implicit none

        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: sp_Cptr

        ! local variables
        type(StrAer), pointer                   :: aer
        type(StrDim), pointer                 :: dimen
        type(StrSpecData), pointer              :: sp
        
        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(sp_Cptr, sp)
                
        call allocate_aer(aer, dimen, sp)

    end subroutine PS_allocate_aer

    ! C-bound wrapper for deallocate_aer
    subroutine PS_deallocate_aer(aer_Cptr) bind(C, name='PS_deallocate_aer')
        implicit none

        type(c_ptr), value, intent(in)          :: aer_Cptr

        ! local variables
        type(StrAer), pointer                   :: aer
        
        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
                
        call deallocate_aer(aer)

    end subroutine PS_deallocate_aer

    ! C-bound wrapper for allocate_aer_prsc
    subroutine PS_allocate_aer_prsc(aer_Cptr, dimen_Cptr, sp_Cptr) bind(C, name='PS_allocate_aer_prsc')
        implicit none

        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: sp_Cptr
        ! local variables
        type(StrAer), pointer                   :: aer
        type(StrDim), pointer                 :: dimen
        type(StrSpecData), pointer              :: sp
        
        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(sp_Cptr, sp)
                
        call allocate_aer_prsc(aer, dimen, sp)

    end subroutine PS_allocate_aer_prsc
   
    ! C-bound wrapper for deallocate_aer_prsc
    subroutine PS_deallocate_aer_prsc(aer_Cptr) bind(C, name='PS_deallocate_aer_prsc')
        implicit none

        type(c_ptr), value, intent(in)          :: aer_Cptr

        ! local variables
        type(StrAer), pointer                   :: aer
        
        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
                
        call deallocate_aer_prsc(aer)

    end subroutine PS_deallocate_aer_prsc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrBound
    !   boundary conditions at top-of-atmosphere and surface (input fluxes, albedo/emissivity
    !   etc)
    ! radiance_core/def_bound.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

    ! C-bound wrapper for allocate_bound
    subroutine PS_allocate_bound(bound_Cptr, dimen_Cptr, sp_Cptr) bind(C, name='PS_allocate_bound')
        implicit none

        type(c_ptr), value, intent(in)          :: bound_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: sp_Cptr

        ! local variables
        type(StrBound), pointer                 :: bound
        type(StrDim), pointer                 :: dimen
        type(StrSpecData), pointer              :: sp
        
        ! convert types
        call C_F_POINTER(bound_Cptr, bound)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(sp_Cptr, sp)
                
        call allocate_bound(bound, dimen, sp)

    end subroutine PS_allocate_bound

    ! C-bound wrapper for deallocate_bound
    subroutine PS_deallocate_bound(bound_Cptr) bind(C, name='PS_deallocate_bound')
        implicit none

        type(c_ptr), value, intent(in)          :: bound_Cptr

        ! local variables
        type(StrBound), pointer                 :: bound
        
        ! convert types
        call C_F_POINTER(bound_Cptr, bound)
                
        call deallocate_bound(bound)

    end subroutine PS_deallocate_bound

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! StrOut
    !   all output variables (fluxes and other diagnostics)
    ! radiance_core/def_out.F90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    ! C-bound wrapper for deallocate_out
    subroutine PS_deallocate_out(radout_Cptr) bind(C, name='PS_deallocate_out')
        implicit none

        type(c_ptr), value, intent(in)          :: radout_Cptr

        ! local variables
        type(StrOut), pointer                   :: radout
        
        ! convert types
        call C_F_POINTER(radout_Cptr, radout)
                
        call deallocate_out(radout)

    end subroutine PS_deallocate_out

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate radiance field
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine PS_radiance_calc( &
        control_Cptr, dimen_Cptr, spectrum_Cptr, atm_Cptr, cld_Cptr, aer_Cptr, bound_Cptr, radout_Cptr &
        ) bind(C, name='PS_radiance_calc')
        
        implicit none

        type(c_ptr), value, intent(in)          :: control_Cptr
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: spectrum_Cptr
        type(c_ptr), value, intent(in)          :: atm_Cptr
        type(c_ptr), value, intent(in)          :: cld_Cptr
        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: bound_Cptr
        type(c_ptr), value, intent(in)          :: radout_Cptr

        ! local variables
        type(StrCtrl), pointer                  :: control
        type(StrDim), pointer                   :: dimen
        type(StrSpecData), pointer              :: spectrum
        type(StrAtm), pointer                   :: atm
        type(StrCld), pointer                   :: cld
        type(StrAer), pointer                   :: aer
        type(StrBound), pointer                 :: bound
        type(StrOut), pointer                   :: radout
        
        ! convert types
        call C_F_POINTER(control_Cptr, control)
        call C_F_POINTER(dimen_Cptr, dimen)
        call C_F_POINTER(spectrum_Cptr, spectrum)
        call C_F_POINTER(atm_Cptr, atm)
        call C_F_POINTER(cld_Cptr, cld)
        call C_F_POINTER(aer_Cptr, aer)
        call C_F_POINTER(bound_Cptr, bound)
        call C_F_POINTER(radout_Cptr, radout)
                
        call radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)

    end subroutine PS_radiance_calc


  end module SOCRATES_C