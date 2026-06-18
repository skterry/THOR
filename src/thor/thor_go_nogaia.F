 
c
c This routine is designed to be the automated reduction
c routine for Sean Terry's GO-17776 program, which will
c be taking precursor observations for the Roman field.
c
c This program performs several steps:
c 1) identify the photmetry files (hst1pass output)
c 2) read in the photometry files
c 3) determine the boundaries
c 4) build internal reference stars from one exposure (NO Gaia catalog file)
c 5) transform 4 lists into the ref frame
c 6) match each exposure to that internal reference (same as GAIA geometry)
c 7) make stacked images in the reference frame
c
c --- thor_go_nogaia.F fork -------------------------------------------------
c Unlike thor_go_gaia.F, STEP4 does NOT read any external Gaia database file.
c Reference stars are taken from exposure NREF (default 1): bright stars with
c m_ln < -11, same cut as STEP6 uses for matching.
c ------------------------------------------------------------------------
 
 
 
 
 
#define _LINUX_  .true.
#define _NSTMAX_ 199999
#define _NMATMX_  99999
 
                                              
c*********************************************
c****                                         
c**** #include "sub_query_GAIA_DR3.f"
c****                                         
c*********************************************
                                              
 
c---------------------------------------------------------------------
c
c 3 million stars in 0605 ; largest number ; 117 files
c
c 1811709771 in the catalog ; agrees with their numbers!
c
c 4*1801*3601 = 25941604/85 = 305195.34   ---> 305201 = the first star
c
c
c   Os:    1811709771
c size:     150386065/85 =
c size:  153995330535/85 = 1811709771
c
c---------------------------------------------------------------------
 
#define _GMAX_  999 999
 
c---------------------------------------------------------------------
c
c this code is designed to be compiled with GFORTRAN:
c
c gfortran query_GAIA_DR3.F -o query_GAIA_DR3.F.e
c
c To execute:
c
c ./query_GAIA_DR3.e <return>
c
c (when you run it with no arguments, the routine will
c  tell you the possible arguments)
c
c---------------------------------------------------------------------
c
c 2019.09.30 - fixed a bug that unselected stars with
c              large errors (no warning); I modified the
c              default limits to hopefully gi-normous
c              levels
c 2018.07.23 - fixed a bug that printed out the wrong
c              "auxiliary" info for stars in cone
c              searches when RV was printed
c 2018.06.12 - for the code mode you now have to enter
c              RA_CEN= DE_CEN= and RAD=
c 2018.06.12 - (historical note from Gaia tools; not used in dex_no_gaia)
c 2018.06.12 - include a mode where we identify stars by ID number
c              this requires a second data base, BASE.GDR3_IDMAP,
c              which has the stars in order, along with the record
c              number where each can be found in the main database.
c 2018.06.08 - redid the catalog; there were a few stars (70) missing for
c              various reasons (bad gz files at the source, some at the
c              edge of dec-strips).  The catalog now has all the stars.
c              Added a few paramters (suggested by Mattia)
c              Added RV.
c              Added the XY mode to push things into an HST-friendly
c              frame.
c 2018.05.22 - implemented the DUMP mode; best to modify
c              the routine query_GAIA_DR2_DUMP subroutine
c              to selectively ouput over the whole sky
c 2018.05.02 - fixed some clipping issues by expanding the blocks searched
c              by one all around
c 2018.05.03 - add flexibility to do a box and not sort
c 2018.05.07 - fixed some minor nmap_ij mapping issues;
c              symptoms were that some squares of the sky didn't
c              get searched;
c              dB and code should be good now
c
 
 
      subroutine query_gaia
      implicit none
 
      real*8  RA_CEN      ! RA center of search
      real*8  DE_CEN      ! Dec center of search
      real*4  RAD         ! degrees
 
      real*8,    dimension(:), allocatable ::    ra_o
      real*8,    dimension(:), allocatable ::    de_o
      real*8,    dimension(:), allocatable ::   era_o
      real*8,    dimension(:), allocatable ::   ede_o
      real*8,    dimension(:), allocatable ::    mg_o
      real*8,    dimension(:), allocatable ::   emg_o
      real*8,    dimension(:), allocatable ::    mb_o
      real*8,    dimension(:), allocatable ::   emb_o
      real*8,    dimension(:), allocatable ::    mr_o
      real*8,    dimension(:), allocatable ::   emr_o
      real*8,    dimension(:), allocatable ::  mura_o
      real*8,    dimension(:), allocatable :: emura_o
      real*8,    dimension(:), allocatable ::  mude_o
      real*8,    dimension(:), allocatable :: emude_o
      real*8,    dimension(:), allocatable ::    pi_o
      real*8,    dimension(:), allocatable ::   epi_o
      integer*8, dimension(:), allocatable ::    id_o
      byte,    dimension(:,:), allocatable ::   b85_o
 
      integer Os, Ou
 
      integer,   dimension(:), allocatable ::  o_o
      real*4,    dimension(:), allocatable ::  v_o
 
      character*80 STRING_PROG
      character*40 STRING_RA, STRING_RAo
      character*40 STRING_DE, STRING_DEo
      character*40 STRING_RAD
      character*20 STRING_DIAG
      character*80 STRING80
      character*80 STRING20
      character*80 FILERDML
 
      logical RA_SEXIG
      logical DE_SEXIG
      logical DE_NEG
 
      real RAh, RAm, RAs
      real DEd, DEm, DEs
 
      integer i
      integer O
      real*8 dra, dde, ddd
      real*8 rd2x, rd2y
      real*8 xy2r, xy2d
 
      integer NLIST, NLISTs
      character*99 LIST_STRING99
      integer NLS99
 
      character*80 FILE_LISTID
      logical DOLISTID
      logical DOLISTRD
      logical DOBOX
      logical DODUMP
 
      logical DOXY
      integer XY_XCEN
      integer XY_YCEN
      real    XY_PSCL
      real    XY_RDAT
      logical SHOWFULL
 
      real*8 dYEARS
      real*8 RA_RDATE, RA_RDATEx
      real*8 DE_RDATE, DE_RDATEx
      real*8 XFRAME, XFRAMEx
      real*8 YFRAME, YFRAMEx
 
      real*8  RA1, RA2, DE1, DE2
      logical DIAG
 
      integer NARG, NARGs
 
      real*8 dclip
 
      integer L, Ls
      integer*8 id_l(_GMAX_), idu
      integer*4 r_l
      integer*4 r_id
 
      real*8   RA_MINU,    RA_MAXU
      real*8   DE_MINU,    DE_MAXU
      real*8  ERA_MINU,   ERA_MAXU
      real*8  EDE_MINU,   EDE_MAXU
      real*8   MG_MINU,    MG_MAXU
      real*8  EMG_MINU,   EMG_MAXU
      real*8   MB_MINU,    MB_MAXU
      real*8  EMB_MINU,   EMB_MAXU
      real*8   MR_MINU,    MR_MAXU
      real*8   EMR_MINU,  EMR_MAXU
      real*8  MURA_MINU, MURA_MAXU
      real*8 EMURA_MINU,EMURA_MAXU
      real*8  MUDE_MINU, MUDE_MAXU
      real*8 EMUDE_MINU,EMUDE_MAXU
      real*8    PI_MINU,   PI_MAXU
      real*8   EPI_MINU,  EPI_MAXU
      real*8  NVIS_MINU, NVIS_MAXU
      real*8  NMAT_MINU, NMAT_MAXU
      real*8  IGOF_MINU, IGOF_MAXU
      real*8  IXSN_MINU, IXSN_MAXU
      real*8  NBAD_MINU, NBAD_MAXU
      real*8  NGUD_MINU, NGUD_MAXU
      real*8  UFLG_MINU, UFLG_MAXU
      real*8  GFLG_MINU, GFLG_MAXU
      real*8  RFLG_MINU, RFLG_MAXU
      real*8    RV_MINU,   RV_MAXU
      real*8   ERV_MINU,  ERV_MAXU
      real*8   BmR_MINU,  BmR_MAXU
      real*8   BmG_MINU,  BmG_MAXU
      real*8   GmR_MINU,  GmR_MAXU
 
      data   RA_MINU,    RA_MAXU / -999999.99d0, 999999.99d0 /
      data   DE_MINU,    DE_MAXU / -999999.99d0, 999999.99d0 /
      data  ERA_MINU,   ERA_MAXU / -999999.99d0, 999999.99d0 /
      data  EDE_MINU,   EDE_MAXU / -999999.99d0, 999999.99d0 /
      data   MG_MINU,    MG_MAXU / -999999.99d0, 999999.99d0 /
      data  EMG_MINU,   EMG_MAXU / -999999.99d0, 999999.99d0 /
      data   MB_MINU,    MB_MAXU / -999999.99d0, 999999.99d0 /
      data  EMB_MINU,   EMB_MAXU / -999999.99d0, 999999.99d0 /
      data   MR_MINU,    MR_MAXU / -999999.99d0, 999999.99d0 /
      data   EMR_MINU,  EMR_MAXU / -999999.99d0, 999999.99d0 /
      data  MURA_MINU, MURA_MAXU / -999999.99d0, 999999.99d0 /
      data EMURA_MINU,EMURA_MAXU / -999999.99d0, 999999.99d0 /
      data  MUDE_MINU, MUDE_MAXU / -999999.99d0, 999999.99d0 /
      data EMUDE_MINU,EMUDE_MAXU / -999999.99d0, 999999.99d0 /
      data    PI_MINU,   PI_MAXU / -999999.99d0, 999999.99d0 /
      data   EPI_MINU,  EPI_MAXU / -999999.99d0, 999999.99d0 /
      data  NVIS_MINU, NVIS_MAXU / -999999.99d0, 999999.99d0 /
      data  NMAT_MINU, NMAT_MAXU / -999999.99d0, 999999.99d0 /
      data  IGOF_MINU, IGOF_MAXU / -999999.99d0, 999999.99d0 /
      data  IXSN_MINU, IXSN_MAXU / -999999.99d0, 999999.99d0 /
      data  NBAD_MINU, NBAD_MAXU / -999999.99d0, 999999.99d0 /
      data  NGUD_MINU, NGUD_MAXU / -999999.99d0, 999999.99d0 /
      data  UFLG_MINU, UFLG_MAXU / -999999.99d0, 999999.99d0 /
      data  GFLG_MINU, GFLG_MAXU / -999999.99d0, 999999.99d0 /
      data  RFLG_MINU, RFLG_MAXU / -999999.99d0, 999999.99d0 /
      data    RV_MINU,   RV_MAXU / -999999.99d0, 999999.99d0 /
      data   ERV_MINU,  ERV_MAXU / -999999.99d0, 999999.99d0 /
      data   BmG_MINU,  BmG_MAXU / -999999.99d0, 999999.99d0 /
      data   BmR_MINU,  BmR_MAXU / -999999.99d0, 999999.99d0 /
      data   GmR_MINU,  GmR_MAXU / -999999.99d0, 999999.99d0 /
 
      common / LIMITS_ /  RA_MINU,    RA_MAXU,
     .                    DE_MINU,    DE_MAXU,
     .                   ERA_MINU,   ERA_MAXU,
     .                   EDE_MINU,   EDE_MAXU,
     .                    MG_MINU,    MG_MAXU,
     .                   EMG_MINU,   EMG_MAXU,
     .                    MB_MINU,    MB_MAXU,
     .                   EMB_MINU,   EMB_MAXU,
     .                    MR_MINU,    MR_MAXU,
     .                    EMR_MINU,  EMR_MAXU,
     .                   MURA_MINU, MURA_MAXU,
     .                  EMURA_MINU,EMURA_MAXU,
     .                   MUDE_MINU, MUDE_MAXU,
     .                  EMUDE_MINU,EMUDE_MAXU,
     .                     PI_MINU,   PI_MAXU,
     .                    EPI_MINU,  EPI_MAXU,
     .                   NVIS_MINU, NVIS_MAXU,
     .                   NMAT_MINU, NMAT_MAXU,
     .                   IGOF_MINU, IGOF_MAXU,
     .                   IXSN_MINU, IXSN_MAXU,
     .                   NBAD_MINU, NBAD_MAXU,
     .                   NGUD_MINU, NGUD_MAXU,
     .                   UFLG_MINU, UFLG_MAXU,
     .                   GFLG_MINU, GFLG_MAXU,
     .                   RFLG_MINU, RFLG_MAXU,
     .                     RV_MINU,   RV_MAXU,
     .                    ERV_MINU,  ERV_MAXU,
     .                    BmG_MINU,  BmG_MAXU,
     .                    BmR_MINU,  BmR_MAXU,
     .                    GmR_MINU,  GmR_MAXU
 
      real*8 date
      real*8 ra_rad
      real*8 de_rad
      real*8 par_ra
      real*8 par_de
 
      real*8 sig_ra, SIGRA_FRAME
      real*8 sig_de, SIGDE_FRAME
 
      integer UNIT
      character*80 FULL_FILE
 
      real ruwe_id
 
 
      FULL_FILE = 'NONE'
 
      ra_rad =  90.00000
      de_rad = -66.56071
 
      ra_rad =   5.65833378
      de_rad = -72.06497192
 
c
c      do i = 0, 1000
c         date = 2000 + i/365.25
c         call calcpi(ra_rad,de_rad,date,par_ra,par_de)
c         write(*,219) i, date, ra_rad, de_rad,
c     .                         par_ra, par_de
c  219    format(1x,i4,1x,f10.4,
c     .          5x,f12.8,1x,f12.8,
c     .          5x,f12.8,1x,f12.8)
c         enddo
c
c     if (.true.) stop
c
 
 
      RA_CEN = 90.00
      DE_CEN =  0.00
      RAD    =  0.20
      STRING_RAD = '0.01'
 
      DOLISTRD  = .false.
      DOLISTID  = .false.
      DOBOX     = .false.
      DODUMP    = .false.
      DOXY      = .false.
      DIAG      = .true.
      DIAG      = .false.
      SHOWFULL  = .false.
 
      allocate(   ra_o(_GMAX_))
      allocate(   de_o(_GMAX_))
      allocate(  era_o(_GMAX_))
      allocate(  ede_o(_GMAX_))
      allocate(   mg_o(_GMAX_))
      allocate(  emg_o(_GMAX_))
      allocate(   mb_o(_GMAX_))
      allocate(  emb_o(_GMAX_))
      allocate(   mr_o(_GMAX_))
      allocate(  emr_o(_GMAX_))
      allocate( mura_o(_GMAX_))
      allocate(emura_o(_GMAX_))
      allocate( mude_o(_GMAX_))
      allocate(emude_o(_GMAX_))
      allocate(   pi_o(_GMAX_))
      allocate(  epi_o(_GMAX_))
      allocate(   id_o(_GMAX_))
      allocate(  b85_o(85,_GMAX_))
 
      allocate(    o_o(_GMAX_))
      allocate(    v_o(_GMAX_))
 
      if (iargc().eq.0) then
         print*,' '
         print*,'There are currently FIVE ways to run:   '
         print*,' '
         print*,'For any, if you provide (on the command line) an '
         print*,'output file name like this:   '
         print*,' '
         print*,'    FULL_FILE=fuller_dr3.out'
         print*,' '
         print*,'the routine will output not just the basic info, '
         print*,'but also the more detailed info, such as the     '
         print*,'correlation info and the astrometric and photometric'
         print*,'quality info.'
         print*,' '
         print*,'-------------------------------------------------'
         print*,'OPTION#1 ---  cone search '
         print*,'-------------------------------------------------'
         print*,' '
         print*,'    3+ args:  RA_CEN=RA DE_CEN=Dec RAD=rad (deg)   '
         print*,'              [DIAG+]  '
         print*,'              [RV+] '
         print*,'              ["XY(XCEN,YCEN,PIX_SCL,RDATE)"]'
         print*,' '
         print*,'                      with '
         print*,' '
         print*,'                      RA  - SEXIG or DECIMAL '
         print*,'                      Dec - SEXIG or DECIMAL '
         print*,'                                             '
         print*,' '
         print*,'examples: '
         print*,' '
         print*,'./query_GAIA_DR3.e "RA_CEN=00:24:05.359"      \'
         print*,'                   "DE_CEN=-72:04:53.20"      \'
         print*,'                   "RAD=1.0"                  \'
         print*,'                   "FULL_FILE=fuller_dr3.out"  '
         print*,' '
         print*,'./query_GAIA_DR2.e "RA_CEN= 24.335"           \'
         print*,'                   "DE_CEN=-72.103            \'
         print*,'                    RAD=0.075                 \'
         print*,'                   "XY(5000,5000,50,2017.3)"   '
         print*,'                                                      '
         print*,'                    (make a frame with the ACS SCALE) '
         print*,'                    (needs to be in quotes to prevent '
         print*,'                                        parser error) '
         print*,'                    (I have not yet double-checked the'
         print*,'                        time-pushing; need to add PI) '
         print*,'                                                      '
         print*,'In the standard-out output, the routine also reports  '
         print*,'the star number (in this local list) and the          '
         print*,'tangent-plane offset in RA, Dec and the total offset  '
         print*,'(all in arcseconds).                                  '
         print*,'                                                      '
         print*,'------------------------------------------------------'
         print*,'OPTION#2 ---                                          '
         print*,'------------------------------------------------------'
         print*,'                                                      '
         print*,'   ./query_GAIA.e DUMP                                '
         print*,'                                                      '
         print*,'This can work like the BOX search did before.  W/o    '
         print*,'delimeters, it just outputs the whole catalog to      '
         print*,'standard out.  It selects on quantities, so you can   '
         print*,'search over a box search of RA/DEC.  You can select   '
         print*,'(say) stars with blue colors or high PM or high       '
         print*,'parallax.  It takes about an hour to search through   '
         print*,'the whole catalog, longer if it has to output a a     '
         print*,'large number of points.                               '
         print*,'                                                      '
         print*,'This outputs the all the available columns, not just  '
         print*,'the "main" ones.                                      '
         print*,'                                                      '
         print*,'------------------------------------------------------'
         print*,'OPTION#3 ---                                          '
         print*,'                                                      '
         print*,'   ./query_GAIA.e LISTRD=list.rd RAD=0.01             '
         print*,'                                                      '
         print*,'This option will take a list of (RA,Dec) coordinates  '
         print*,'in a file (in either decimal or sexigesimal format)   '
         print*,'and it will do a cone search around each one,         '
         print*,'outputting the basic info to the screen.              '
         print*,'                                                      '
         print*,'------------------------------------------------------'
         print*,'                                                      '
         print*,'THE REST OF THE OPTIONS ARE NOT WORKING YET; I CAN    '
         print*,'GET TO THEM IF THERE IS AN INTEREST.  JUST  LET ME    '
         print*,'KNOW.                                                 '
         print*,'                                                      '
         print*,'OPTION#4 --- by list ID                               '
         print*,'                                                      '
         print*,'------------------------------------------------------'
         print*,' '
         print*,'List of quantities you can set delimeters on: '
         print*,' '
         print*,'        RA_MINU   RA_MAXU   ERA_MINU   ERA_MAXU'
         print*,'        DE_MINU   DE_MAXU   EDE_MINU   EDE_MAXU'
         print*,'        MG_MINU   MG_MAXU   EMG_MINU   EMG_MAXU'
         print*,'        MB_MINU   MB_MAXU   EMB_MINU   EMB_MAXU'
         print*,'        MR_MINU   MR_MAXU   EMR_MINU   EMR_MAXU'
         print*,'      MURA_MINU MURA_MAXU EMURA_MINU EMURA_MAXU'
         print*,'      MUDE_MINU MUDE_MAXU EMUDE_MINU EMUDE_MAXU'
         print*,'        PI_MINU   PI_MAXU   EPI_MINU   EPI_MAXU'
         print*,' '
         print*,'      And colors:  '
         print*,'       BmG_MINU  BmG_MAXU   BmR_MINU   BmR_MAXU'
         print*,'                            GmR_MINU   GmR_MAXU'
         print*,' '
         print*,' '
         print*,' '
         stop
         endif
 
      NARGs = iargc()
 
      call getarg(0,STRING_PROG)
      write(*,'(''# ARG'',i2.2,'': '',80a)') 0, STRING_PROG
      do NARG = 1, NARGs
         call getarg(NARG,STRING80)
         write(*,'(''# ARG'',i2.2,'': '',80a)') NARG, STRING80
         if (STRING80(1:3).eq.'XY(') then
             DOXY = .true.
             do i = 4, 80
                if (STRING80(i:i).eq.')') STRING80(i:i) = ' '
                if (STRING80(i:i).eq.',') STRING80(i:i) = ' '
                enddo
             read(STRING80(4:80),*) XY_XCEN, XY_YCEN, XY_PSCL, XY_RDAT
             if (XY_RDAT.lt.1950.or.XY_RDAT.gt.2050)
     .           stop 'XY_RDAT should be in the form of 2012.33'
             goto 1
             endif
         if (STRING80(1:10).eq.'FULL_FILE=') then
             FULL_FILE = STRING80(11:80)
             SHOWFULL = .true.
             open(50,file=FULL_FILE,status='unknown')
             goto 1
             endif
         if (STRING80(1:5).eq.'DIAG+') then
             DIAG = .true.
             goto 1
             endif
         if (STRING80(1:5).eq.'DIAG-') then
             DIAG = .false.
             goto 1
             endif
         if (STRING80(1:8).eq.'RA_MINU=') then
             read(STRING80(9:80),*) RA_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'RA_MAXU=') then
             read(STRING80(9:80),*) RA_MAXU
             goto 1
             endif
         if (STRING80(1:8).eq.'DE_MINU=') then
             read(STRING80(9:80),*) DE_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'DE_MAXU=') then
             read(STRING80(9:80),*) DE_MAXU
             goto 1
             endif
         if (STRING80(1:8).eq.'MG_MINU=') then
             read(STRING80(9:80),*) MG_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'MG_MAXU=') then
             read(STRING80(9:80),*) MG_MAXU
             goto 1
             endif
         if (STRING80(1:8).eq.'MR_MINU=') then
             read(STRING80(9:80),*) MR_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'MR_MAXU=') then
             read(STRING80(9:80),*) MR_MAXU
             goto 1
             endif
         if (STRING80(1:8).eq.'MB_MINU=') then
             read(STRING80(9:80),*) MB_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'MB_MAXU=') then
             read(STRING80(9:80),*) MB_MAXU
             goto 1
             endif
         if (STRING80(1:8).eq.'PI_MINU=') then
             read(STRING80(9:80),*) PI_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'PI_MAXU=') then
             read(STRING80(9:80),*) PI_MAXU
             goto 1
             endif
         if (STRING80(1:8).eq.'RV_MINU=') then
             read(STRING80(9:80),*) RV_MINU
             goto 1
             endif
         if (STRING80(1:8).eq.'RV_MAXU=') then
             read(STRING80(9:80),*) RV_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'ERA_MAXU=') then
             read(STRING80(10:80),*) ERA_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'ERA_MINU=') then
             read(STRING80(10:80),*) ERA_MINU
             goto 1
             endif
         if (STRING80(1:9).eq.'EDE_MAXU=') then
             read(STRING80(10:80),*) EDE_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'EDE_MINU=') then
             read(STRING80(10:80),*) EDE_MINU
             goto 1
             endif
         if (STRING80(1:9).eq.'EMG_MAXU=') then
             read(STRING80(10:80),*) EMG_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'EMG_MINU=') then
             read(STRING80(10:80),*) EMG_MINU
             goto 1
             endif
         if (STRING80(1:9).eq.'EMB_MAXU=') then
             read(STRING80(10:80),*) EMB_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'EMB_MINU=') then
             read(STRING80(10:80),*) EMB_MINU
             goto 1
             endif
         if (STRING80(1:9).eq.'EMR_MAXU=') then
             read(STRING80(10:80),*) EMR_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'EMR_MINU=') then
             read(STRING80(10:80),*) EMR_MINU
             goto 1
             endif
         if (STRING80(1:9).eq.'EPI_MAXU=') then
             read(STRING80(10:80),*) EPI_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'EPI_MINU=') then
             read(STRING80(10:80),*) EPI_MINU
             goto 1
             endif
         if (STRING80(1:9).eq.'ERV_MAXU=') then
             read(STRING80(10:80),*) ERV_MAXU
             goto 1
             endif
         if (STRING80(1:9).eq.'ERV_MINU=') then
             read(STRING80(10:80),*) ERV_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'MURA_MINU=') then
             read(STRING80(11:80),*) MURA_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'MURA_MAXU=') then
             read(STRING80(11:80),*) MURA_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'MUDE_MINU=') then
             read(STRING80(11:80),*) MUDE_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'MUDE_MAXU=') then
             read(STRING80(11:80),*) MUDE_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'NVIS_MINU=') then
             read(STRING80(11:80),*) NVIS_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'NVIS_MAXU=') then
             read(STRING80(11:80),*) NVIS_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'NMAT_MINU=') then
             read(STRING80(11:80),*) NMAT_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'NMAT_MAXU=') then
             read(STRING80(11:80),*) NMAT_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'IGOF_MINU=') then
             read(STRING80(11:80),*) IGOF_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'IGOF_MAXU=') then
             read(STRING80(11:80),*) IGOF_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'IXSN_MINU=') then
             read(STRING80(11:80),*) IXSN_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'IXSN_MAXU=') then
             read(STRING80(11:80),*) IXSN_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'NBAD_MINU=') then
             read(STRING80(11:80),*) NBAD_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'NBAD_MAXU=') then
             read(STRING80(11:80),*) NBAD_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'NGUD_MINU=') then
             read(STRING80(11:80),*) NGUD_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'NGUD_MAXU=') then
             read(STRING80(11:80),*) NGUD_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'UFLG_MINU=') then
             read(STRING80(11:80),*) UFLG_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'UFLG_MAXU=') then
             read(STRING80(11:80),*) UFLG_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'GFLG_MINU=') then
             read(STRING80(11:80),*) GFLG_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'GFLG_MAXU=') then
             read(STRING80(11:80),*) GFLG_MAXU
             goto 1
             endif
         if (STRING80(1:10).eq.'RFLG_MINU=') then
             read(STRING80(11:80),*) RFLG_MINU
             goto 1
             endif
         if (STRING80(1:10).eq.'RFLG_MAXU=') then
             read(STRING80(11:80),*) RFLG_MAXU
             goto 1
             endif
         if (STRING80(1:11).eq.'EMURA_MINU=') then
             read(STRING80(12:80),*) EMURA_MINU
             goto 1
             endif
         if (STRING80(1:11).eq.'EMURA_MAXU=') then
             read(STRING80(12:80),*) EMURA_MAXU
             goto 1
             endif
         if (STRING80(1:11).eq.'EMUDE_MINU=') then
             read(STRING80(12:80),*) EMUDE_MINU
             goto 1
             endif
         if (STRING80(1:11).eq.'EMUDE_MAXU=') then
             read(STRING80(12:80),*) EMUDE_MAXU
             goto 1
             endif
         if (STRING80(1:7).eq.'LISTRD=') goto 1
         if (STRING80(1:7).eq.'LISTID=') goto 1
         if (STRING80(1:4).eq.'BOX=')    goto 1
         if (STRING80(1:4).eq.'DUMP')    goto 1
         if (STRING80(1:4).eq.'DUMP')    goto 1
         if (STRING80(1:7).eq.'RA_CEN=') goto 1
         if (STRING80(1:7).eq.'DE_CEN=') goto 1
         if (STRING80(1:4).eq.'RAD=') then
             STRING_RAD = STRING80(5:40)
             goto 1
             endif
         print*,'PROBLEM PARSING: ',STRING80
         stop
    1    continue
         enddo
 
 
      call getarg(0,STRING_PROG)
      call getarg(1,STRING80)
 
 
      if (STRING80(1:7).eq.'LISTID=') then
          stop 'NOT READY FOR LISTID YET'
          endif
 
 
      if (STRING80(1:7).eq.'LISTRD=') then
          if (DOXY) stop 'Cannot do LIST and DOXY.'
          DOLISTRD = .true.
          read(STRING80(8:80),*) FILERDML
          open(21,file=FILERDML,status='old')
          endif
 
 
      if (STRING80(1:4).eq.'DUMP') then
          call sub_query_GAIA_DR3_DUMP(RA_MINU,RA_MAXU,
     .                                 DE_MINU,DE_MAXU)
          stop
          endif
 
 
      if (.not.DOLISTRD) then
          call getarg(1,STRING_RAo)
          call getarg(2,STRING_DEo)
          call getarg(3,STRING_RAD)
          if (STRING_RAo(1:7).ne.'RA_CEN='.or.
     .        STRING_DEo(1:7).ne.'DE_CEN='.or.
     .        STRING_RAD(1:4).ne.'RAD=') then
              print*,'                                   '
              print*,'If you are not running this in     '
              print*,'LIST mode, BOX mode, or DUMP mode, '
              print*,'then you must want CONE mode.      '
              print*,'This requires the first three args '
              print*,'to be RA_CEN=... DE_CEN=... RAD=.. '
              print*,'                                   '
              stop
              endif
          STRING_RAo = STRING_RAo(8:40)
          STRING_DEo = STRING_DEo(8:40)
          STRING_RAD = STRING_RAD(5:40)
          call show_delimeters
          endif
 
      NLIST = 0
 333  continue
      if (DOLISTRD) then
          read(21,'(a99)',end=334) LIST_STRING99
          read(LIST_STRING99,*) STRING_RAo, STRING_DEo
          NLIST = NLIST+1
          NLS99 = 1
          do i = 1, 99
             if (LIST_STRING99(i:i).ne.' ') NLS99 = i
             enddo
          endif
 
      STRING_RA = STRING_RAo
      STRING_DE = STRING_DEo
 
      RA_SEXIG = .false.
      do i = 1, 20
         if (STRING_RA(i:i).eq.':') then
             RA_SEXIG = .true.
             STRING_RA(i:i) = ' '
             endif
         enddo
      if (RA_SEXIG) then
          read(STRING_RA,*) RAh, RAm, RAs
          RA_CEN = (RAh + RAm/60.0d0 + RAs/60.0d0/60.0d0)/24.0d0*360.0d0
         else
          read(STRING_RA,*) RA_CEN
         endif
 
      DE_NEG   = .false.
      DE_SEXIG = .false.
      do i = 1, 20
         if (STRING_DE(i:i).eq.'-') then
             DE_NEG = .true.
             STRING_DE(i:i) = ' '
             endif
         if (STRING_DE(i:i).eq.':') then
             DE_SEXIG = .true.
             STRING_DE(i:i) = ' '
             endif
         enddo
      if (DE_SEXIG) then
          read(STRING_DE,*) DEd, DEm, DEs
          DE_CEN = (DEd + DEm/60.0d0 + DEs/60.0d0/60.0d0)
         else
          read(STRING_DE,*) DE_CEN
         endif
      if (DE_NEG) DE_CEN = -DE_CEN
 
      read(STRING_RAD,*) RAD
 
c     print*,'enter sub_query...'
      call sub_query_GAIA_DR3(RA_CEN,DE_CEN,RAD,
     .                        ra_o,  de_o, era_o, ede_o,
     .                        mg_o, emg_o,
     .                        mb_o, emb_o,
     .                        mr_o, emr_o,
     .                      mura_o, emura_o,
     .                      mude_o, emude_o,
     .                        pi_o, epi_o,
     .                        id_o,
     .                       b85_o,
     .                        Os,DIAG)
c     print*,'exit sub_query...',Os
 
 
      if (DIAG) then
         open(88,file='query_GAIA_DR2.diag.05.Os',status='unknown')
         do O = 1, Os
            write(88,188) ra_o(O), de_o(O),
     .                    int(era_o(O)), int(ede_o(O)),
     .                    mg_o(O), int(emg_o(O)),
     .                    mb_o(O), int(emb_o(O)),
     .                    mr_o(O), int(emr_o(O)),
     .                    mura_o(O), emura_o(O),
     .                    mude_o(O), emude_o(O),
     .                    pi_o(O),   epi_o(O),
     .                    id_o(O)
  188       format(f14.10,1x,f14.10,1x,i5,1x,i5,1x,
     .             f6.3,1x,i3,1x,
     .             f6.3,1x,i3,1x,
     .             f6.3,1x,i3,1x,
     .             f8.3,1x,f5.3,1x,
     .             f8.3,1x,f5.3,1x,
     .             f8.3,1x,f5.3,1x,
     .             i19,1x,
     .             i3,1x,i3,1x,
     .             f6.3,1x,
     .             f6.3,1x,
     .             i4,1x,i4,1x,
     .             3i1,1x,
     .             i6,1x,
     .             f8.3,1x,f6.3)
            enddo
         close(88)
         endif
 
 
 
      if (Os.gt._GMAX_) then
          print*,' '
          print*,' THE SEARCH RETURNED TOO MANY SOURCES.'
          print*,' IT IS HARD-CODED TO RETURN AT MOST '
          print*,' _GMAX_ = ',_GMAX_
          print*,' '
          print*,' BUT THE SEARCH RESULTED IN'
          print*,'    Ns  = ',Os
          print*,' '
          print*,' YOU NEED TO INCREASE FINDING LIMIT (_GMAX_).'
          print*,' '
          stop
          endif
 
 
c----------------------------------------------
c
c output the header...
c
      if (.not.DOLISTRD.or.NLIST.eq.1) then
         write(*,'(''#  '')')
         write(*,'(''# ROUTINE: '',a80)') STRING_PROG
         if (.not.DOLISTRD) then
            write(*,'(''#          '',a20)') STRING_RAo
            write(*,'(''#          '',a20)') STRING_DEo
            write(*,'(''#          '',a20)') STRING_RAD
            write(*,'(''#  '')')
            write(*,'(''# RA_CEN: '',f18.8)') RA_CEN
            write(*,'(''# DE_CEN: '',f18.8)') DE_CEN
            write(*,'(''#    RAD: '',f18.8)') RAD
            write(*,'(''#  '')')
            write(*,'(''#      #: '',i12)')   Os
            write(*,'(''#  '')')
            if (DOXY) then
                write(*,'(''#  '')')
                write(*,'(''# PUT OUTPUT INTO A HST-TYPE FRAME:  '')')
                write(*,'(''#        XY_XCEN: '',i5  )') XY_XCEN
                write(*,'(''#        XY_YCEN: '',i5  )') XY_YCEN
                write(*,'(''#        XY_PSCL: '',f8.2)') XY_PSCL
                write(*,'(''#        XY_RDAT: '',f8.2)') XY_RDAT
                write(*,'(''#  '')')
                call sub_calcpi(RA_CEN,DE_CEN,
     .                          XY_RDAT*1.0d0,
     .                          par_ra,  par_de)
                write(*,'(''#  '')')
                write(*,'(''#  KAPLAN_PIFACs:  '',3f12.6)') par_ra,
     .                                                      par_de,
     .                                                      XY_RDAT
                call find_pifac_USNO(RA_CEN,DE_CEN,
     .                               XY_RDAT*1.0d0,
     .                               par_ra,par_de)
                write(*,'(''#    USNO_PIFACs:  '',3f12.6)') par_ra,
     .                                                      par_de,
     .                                                      XY_RDAT
                write(*,'(''#  '')')
                endif
            endif
         if (DOLISTRD) then
            write(*,'(''#    LIST: '',a80)') STRING80
            endif
         write(*,'(''#  '')')
         write(*,'(''# If you use this for your science, see the '')')
         write(*,'(''# way to credit the catalog at the bottom of'')')
         write(*,'(''# the file.'')')
         write(*,'(''#  '')')
         endif
 
      do O = 1, Os
         o_o(O) = O
         dra = rd2x(ra_o(O),de_o(O),RA_CEN,DE_CEN)
         dde = rd2y(ra_o(O),de_o(O),RA_CEN,DE_CEN)
         ddd = sqrt(dra**2+dde**2)
         v_o(O) = ddd
         enddo
 
      call QSORTI (o_o,Os,v_o)
 
      if (DOXY) then
          write(*,'(''#'')')
          write(*,316)
          write(*,'(''#'')')
          goto 444
          endif
 
      if (DOLISTRD) then
          write(*,117)
          write(*,116)
          write(*,115)
          write(*,117)
          goto 444
          endif
 
      write(*,117)
      write(*,116)
      write(*,115)
      write(*,117)
 
 444  continue
 
c     if (DOLISTRD) Os = 1
 
      do O = 1, Os
         Ou = o_o(O)
         dra = rd2x(ra_o(Ou),de_o(Ou),RA_CEN,DE_CEN)
         dde = rd2y(ra_o(Ou),de_o(Ou),RA_CEN,DE_CEN)
         ddd = sqrt(dra**2+dde**2)
         if (mg_o(Ou).lt.1) emg_o(Ou) = 0.
         if (mb_o(Ou).lt.1) emb_o(Ou) = 0.
         if (mr_o(Ou).lt.1) emr_o(Ou) = 0.
         if (.not.( ddd    .lt.99.9)) ddd = 99.9
         if (DOXY) then
 
c            write(STRING20,'(i19)') id_o(Ou)
c            if (STRING20.eq.' 5853498713160606720') then      ! PROX CEN...
c                   pi_o(Ou) =   768.5003653333918
c                  epi_o(Ou) =     0.20302126916720853 ! MAS
c                 mura_o(Ou) = -3781.3059409680022
c                emura_o(Ou) =     0.08518122633177506
c                 mude_o(Ou) =   769.7656533486586
c                emude_o(Ou) =     0.18895487215250334
c                endif
 
             dYEARS = XY_RDAT - 2016.0d0
 
             call find_pifac_USNO(ra_o(Ou),de_o(Ou),
     .                            XY_RDAT*1.0d0,
     .                            par_ra,par_de)
 
             dra = dYEARS*mura_o(Ou)/1000./3600.
     .           +   pi_o(Ou)*par_ra/1000./3600.
             dde = dYEARS*mude_o(Ou)/1000./3600.
     .           +   pi_o(Ou)*par_de/1000./3600.
             RA_RDATE = xy2r(dra,dde,ra_o(Ou),de_o(Ou))
             DE_RDATE = xy2d(dra,dde,ra_o(Ou),de_o(Ou))
 
             dra = rd2x(RA_RDATE,DE_RDATE,RA_CEN,DE_CEN)
             dde = rd2y(RA_RDATE,DE_RDATE,RA_CEN,DE_CEN)
             XFRAME = XY_XCEN - dra*3600.0d0*1000.0d0/XY_PSCL
             YFRAME = XY_YCEN + dde*3600.0d0*1000.0d0/XY_PSCL
 
             RA_RDATEx = ra_o(Ou) + dYEARS*mura_o(Ou)/1000./3600.
     .                             /cos(de_o(Ou)/180.0*3.1415927)
             DE_RDATEx = de_o(Ou) + dYEARS*mude_o(Ou)/1000./3600.
             dra = rd2x(RA_RDATEx,DE_RDATEx,RA_CEN,DE_CEN)
             dde = rd2y(RA_RDATEx,DE_RDATEx,RA_CEN,DE_CEN)
             XFRAMEx = XY_XCEN - dra*3600.0d0*1000.0d0/XY_PSCL
             YFRAMEx = XY_YCEN + dde*3600.0d0*1000.0d0/XY_PSCL
 
c            print*,'RA_RDATE: ',RA_RDATE,RA_RDATEx,RA_RDATE-RA_RDATEx
c            print*,'DE_RDATE: ',DE_RDATE,DE_RDATEx,DE_RDATE-DE_RDATEx
 
             if (XFRAME.lt.-9999.0) XFRAME = -9999.0
             if (YFRAME.lt.-9999.0) YFRAME = -9999.0
             if (XFRAME.gt.99999.0) XFRAME = 99999.0
             if (YFRAME.gt.99999.0) YFRAME = 99999.0
             sig_ra = sqrt((era_o(Ou)/1000.)**2     ! sigma in RA in mas
     .                    +(dYEARS*emura_o(Ou))**2
     .                    +(par_ra*epi_o(Ou))**2)
             sig_de = sqrt((ede_o(Ou)/1000.)**2     ! sigma in RA in mas
     .                    +(dYEARS*emude_o(Ou))**2
     .                    +(par_de*epi_o(Ou))**2)
             SIGRA_FRAME = sig_ra/XY_PSCL
             SIGDE_FRAME = sig_de/XY_PSCL
 
c            print*,'    pi_o(Ou): ',  pi_o(Ou),   epi_o(Ou),
c    .                                 par_ra, par_de
c            print*,'  mura_o(Ou): ',mura_o(Ou), emura_o(Ou)
c            print*,'  mude_o(Ou): ',mude_o(Ou), emude_o(Ou)
c            print*,'                    '
c            print*,'          era_o(Ou): ',era_o(Ou)/1000.
c            print*,' dYEARS*emura_o(Ou): ',dYEARS*emura_o(Ou)
c            print*,'   par_ra*epi_o(Ou): ',par_ra*epi_o(Ou)
 
             if (abs(dYEARS).le.0.1.or.
     .                (abs(mura_o(Ou)).gt.0.001.and.
     .                 abs(mude_o(Ou)).gt.0.001))
     .            write(*,319) XFRAME, YFRAME, mg_o(Ou),
     .                         SIGRA_FRAME,
     .                         SIGDE_FRAME,
     .                         XFRAMEx, YFRAMEx,
     .                         XY_RDAT,
     .                         RA_RDATE, DE_RDATE,
     .                         ra_o(Ou),        de_o(Ou),
     .                         era_o(Ou)/1000., ede_o(Ou)/1000.,
     .                         mg_o(Ou), emg_o(Ou)/1000.,
     .                         mb_o(Ou), emb_o(Ou)/1000.,
     .                         mr_o(Ou), emr_o(Ou)/1000.,
     .                         mura_o(Ou), emura_o(Ou),
     .                         mude_o(Ou), emude_o(Ou),
     .                         pi_o(Ou),   epi_o(Ou),
     .                         id_o(Ou),  ruwe_id(id_o(Ou))
             goto 888
             endif
         if (DOLISTRD) then
             write(*,119)  ra_o(Ou),        de_o(Ou),
     .                    era_o(Ou)/1000., ede_o(Ou)/1000.,
     .                    mg_o(Ou), emg_o(Ou)/1000.,
     .                    mb_o(Ou), emb_o(Ou)/1000.,
     .                    mr_o(Ou), emr_o(Ou)/1000.,
     .                    mura_o(Ou), emura_o(Ou),
     .                    mude_o(Ou), emude_o(Ou),
     .                    pi_o(Ou),   epi_o(Ou),
     .                    id_o(Ou),
     .                    dclip(dra*3600,-9999.9d0,9999.9d0),
     .                    dclip(dde*3600,-9999.9d0,9999.9d0),
     .                    dclip(ddd*3600,    0.0d0,9999.9d0),
     .                    NLIST
c    .                    nvis_o(Ou), nmat_o(Ou),
c    .                    igof_o(Ou), ixsn_o(Ou),
c    .                    nbad_o(Ou), ngud_o(Ou),
c    .                    uflg_o(Ou), gflg_o(Ou), rflg_o(Ou),
c    .                    rv_o(Ou), erv_o(Ou),
c    .                    RA_CEN, DE_CEN,
c    .                    LIST_STRING99(1:NLS99)
             goto 888
             endif
         write(*,119)  ra_o(Ou),        de_o(Ou),
     .                era_o(Ou)/1000., ede_o(Ou)/1000.,
     .                 mg_o(Ou), emg_o(Ou)/1000.,
     .                 mb_o(Ou), emb_o(Ou)/1000.,
     .                 mr_o(Ou), emr_o(Ou)/1000.,
     .                mura_o(Ou), emura_o(Ou),
     .                mude_o(Ou), emude_o(Ou),
     .                pi_o(Ou),   epi_o(Ou),
     .                id_o(ou),
     .               dclip(dra*3600,-9999.9d0,9999.9d0),
     .               dclip(dde*3600,-9999.9d0,9999.9d0),
     .               dclip(ddd*3600,    0.0d0,9999.9d0),
     .               O, ruwe_id(id_o(Ou))
  888    continue
         if (SHOWFULL) then
             call showfull_b85(b85_o(1,Ou),O.eq.Os,50)
             endif
         enddo
 
      if (SHOWFULL) close(50)
 
  119 format(f14.10,1x,f14.10,1x,
     .       f7.3,1x,f7.3,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f8.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       i19,1x,
     .       f8.2,1x,f8.2,1x,f8.2,1x,i5,f6.2)
  115 format('#    (deg)    ',1x,'      (deg)   ',1x,
     .       '  (mas)',1x,'  (mas)',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '   (mas)',1x,'(mas)',1x,
     .       '                  ',1x,
     .       '   (mas)',1x,'   (mas) ',1x,'  (mas)',1x,
     .       '   (#)',1x,'(qof) ')
  116 format('#RIGHT ASCENSN',1x,'   DECLINATION',1x,
     .       '   eRA ',1x,'   eDEC',1x,
     .       '  GMAG',1x,' eGMAG ',1x,
     .       '  BMAG',1x,' eBMAG ',1x,
     .       '  RMAG',1x,' eRMAG ',1x,
     .       '    muRA',1x,'emuRA',1x,
     .       '   muDEC',1x,'emuDE',1x,
     .       '     PI ',1x,'  ePI',1x,
     .       ' GAIA_DR2_ID_NUMBER ',1x,
     .       '  dRA   ',1x,'  dDE   ',1x,'  dANG ',1x,
     .       '   O',1x,'  RUWE ')
  117 format('#.............',1x,'..............',1x,
     .       '.......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '........',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '...................',1x,
     .       '........',1x,'........',1x,'........',1x,
     .       '.....',1x,'......')
  415 format('#    (deg)    ',1x,'      (deg)   ',1x,
     .       '  (mas)',1x,'  (mas)',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '   (mas)',1x,'(mas)',1x,
     .       '                    ',1x,
     .       '   ',1x,'   ',1x,
     .       '      ',1x,
     .       '      ',1x,
     .       '    ',1x,'    ',1x,
     .       '      ',1x,
     .       '   (km/s)',1x,'(km/s)',1x,
     .       '        ')
  416 format('#RIGHT ASCENSN',1x,'   DECLINATION',1x,
     .       '   eRA ',1x,'   eDEC',1x,
     .       '  GMAG',1x,' eGMAG ',1x,
     .       '  BMAG',1x,' eBMAG ',1x,
     .       '  RMAG',1x,' eRMAG ',1x,
     .       '    muRA',1x,'emuRA',1x,
     .       '   muDEC',1x,'emuDE',1x,
     .       '     PI ',1x,'  ePI',1x,
     .       ' GAIA_DR2_ID_NUMBER ')
c    .       'NVI',1x,'NMT',1x,
c    .       '  IGOF',1x,
c    .       '  IXSN',1x,
c    .       'NBAD',1x,'NGUD',1x,
c    .       ' U G R',1x,
c    .       '      RV',1x,'   eRV',1x,
c    .       '   ORDNO')
c    .       ' ')
  417 format('#.............',1x,'..............',1x,
     .       '.......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '........',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '....................')
c    .       '...',1x,'...',1x,
c    .       '......',1x,
c    .       '......',1x,
c    .       '....',1x,'....',1x,
c    .       ' . . .',1x,
c    .       '........',1x,
c    .       '........',1x,'......')
  319 format(f10.4,1x,f10.4,1x,f7.3,1x,
     .       f10.4,1x,f10.4,1x,
     .       f10.4,1x,f10.4,1x,
     .       f9.4,1x,
     .       f14.10,1x,f14.10,1x,
     .       f14.10,1x,f14.10,1x,
     .       f7.3,1x,f7.3,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f9.3,1x,f5.3,1x,
     .       f9.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       i19,1x,f6.2)
  316 format('#   XFRAME',1x,'    YFRAME',1x,'   GMAG',1x,
     .       '   eXFRAME',1x,'   eYFRAME',1x,
     .       ' XFRAMExPI',1x,' YFRAMExPI',1x,
     .       '    RDATE' ,1x,
     .       '     RA(RDATE)',1x,'    DEC(RDATE)',1x,
     .       '    RA(2016.0)',1x,'   DEC(2016.0)',1x,
     .       '   eRA ',1x,'  eDEC',1x,
     .       '   GMAG',1x,' eGMAG',1x,
     .       '   BMAG',1x,' eBMAG',1x,
     .       '   RMAG',1x,' eRMAG',1x,
     .       '     muRA',1x,' emuRA',
     .       '    muDEC',1x,' emuDE',
     .       '     PI ',1x,'  ePI',2x,
     .       'GAIA_DR2_ID_NUMBER ',
     .       '   RUWE')
  135 format('#    (deg)    ',1x,'      (deg)   ',1x,
     .       '  (mas)',1x,'  (mas)',1x,
     .       ' (mag)',1x,'(mag) ',1x,
     .       ' (mag)',1x,'(mag) ',1x,
     .       ' (mag)',1x,'(mag) ',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '   (mas)',1x,'(mas)',3x,
     .       '                    ')
 
 
      if (DOLISTRD) goto 333
  334 close(21)
 
      if (DOXY) then
         write(*,'(''#'')')
         write(*,316)
         write(*,'(''#'')')
         goto 555
         endif
 
      if (DOLISTRD) then
          write(*,117)
          write(*,116)
          write(*,115)
          write(*,117)
          goto 555
          endif
 
       write(*,117)
       write(*,116)
       write(*,115)
       write(*,117)
 
 555  continue
 
      write(*,'(''# '')' )
      write(*,'(''# -----------------------------------------------'')')
      write(*,'(''# '')' )
      write(*,'(''#    RA -- Right Ascension, degrees      '')')
      write(*,'(''#    DE -- Declination, degrees          '')')
      write(*,'(''#   eRA -- mas; quoted error in RA       '')')
      write(*,'(''#   eDE -- mas; quoted error in RA       '')')
      write(*,'(''#  GMAG -- GAIAs green magnitude         '')')
      write(*,'(''# eGMAG -- error in GMAG (dflx/flx)      '')')
      write(*,'(''#  BMAG -- GAIAs blue magnitude          '')')
      write(*,'(''# eBMAG -- error in BMAG (dflx/flx)      '')')
      write(*,'(''#  RMAG -- GAIAs red magnitude           '')')
      write(*,'(''# eRMAG -- error in BMAG (dflx/flx)      '')')
      write(*,'(''#  muRA -- proper motion in RA (mas/yr)  '')')
      write(*,'(''#          (no cos(DEC) needed)          '')')
      write(*,'(''# emuRA -- err in RA proper motion       '')')
      write(*,'(''#  muDE -- proper motion in DEC (mas/yr) '')')
      write(*,'(''# emuDE -- err in DEC proper motion      '')')
      write(*,'(''#    PI -- parallax (mas)                '')')
      write(*,'(''#   ePI -- err in PI (mas)               '')')
      write(*,'(''#     # -- GAIAs DR2 ID NUMBER (int*8)   '')')
      write(*,'(''#   NVI -- visibility periods used       '')')
      write(*,'(''#   NMT -- num of astrom-matchd obsns    '')')
      write(*,'(''#  IGOF -- astrometric goodness-of-fit   '')')
      write(*,'(''#  IXSN -- astrometric excess noise      '')')
      write(*,'(''# NGOOD -- number of good astrom obsns   '')')
      write(*,'(''#  NBAD -- number of bad astrom obsns    '')')
      write(*,'(''#     U -- flag constructed by Mattia:   '')')
      write(*,'(''#          uflag                         '')')
      write(*,'(''#          = sqrt(chi2/(ngood-5))        '')')
      write(*,'(''#             lt.(1.2*max(1.,            '')')
      write(*,'(''#               exp(-0.2*(gmag-19.5))))  '')')
      write(*,'(''#     G -- flag constructed by Mattia:   '')')
      write(*,'(''#           eval = phot_bp_rp_excess_fac '')')
      write(*,'(''#          gflag = (eval.gt.aux1).and.   '')')
      write(*,'(''#                  (eval.lt.aux2)        '')')
      write(*,'(''#     R -- radial veloicty available?    '')')
      write(*,'(''#    RV -- radial velocity (km/s)        '')')
      write(*,'(''#   eRV -- error in RV (km/s)            '')')
      write(*,'(''#                                                '')')
      write(*,'(''# -----------------------------------------------'')')
      write(*,'(''#                                                '')')
      write(*,'(''# This work has made use of data from the        '')')
      write(*,'(''# European Space Agency (ESA) mission {\it Gaia} '')')
      write(*,'(''# (\url{https://www.cosmos.esa.int/gaia}),       '')')
      write(*,'(''# processed by the {\it Gaia} Data Processing and'')')
      write(*,'(''# Analysis Consortium (DPAC,                     '')')
      write(*,'(''# \url{https://www.cosmos.esa.int/web/gaia/dpac/ '')')
      write(*,'(''#                                 consortium}).  '')')
      write(*,'(''# Funding for the DPAC has been provided by      '')')
      write(*,'(''# national institutions, in particular the       '')')
      write(*,'(''# institutions participating in the {\it Gaia}   '')')
      write(*,'(''# Multilateral Agreement.                        '')')
      write(*,'(''#                                                '')')
      write(*,'(''#                                                '')')
      write(*,'(''# If you have used Gaia DR2 data in your,        '')')
      write(*,'(''# research, please cite both the Gaia mission    '')')
      write(*,'(''# paper and the Gaia DR2 release paper:          '')')
      write(*,'(''#                                                '')')
      write(*,'(''# Gaia Collaboration et al. (2016)               '')')
      write(*,'(''#      "The Gaia mission"  A&A 595, pp. A1.      '')')
      write(*,'(''#           Description of the Gaia mission      '')')
      write(*,'(''#            (spacecraft, instruments, survey and'')')
      write(*,'(''#             meast principles, and operations)  '')')
      write(*,'(''#                                                '')')
      write(*,'(''#  and                                           '')')
      write(*,'(''#                                                '')')
      write(*,'(''# Gaia Collaboration et al. (2018b)              '')')
      write(*,'(''#      "Gaia Data Release 2: Summary of the      '')')
      write(*,'(''#       contents and survey properties.          '')')
      write(*,'(''#       A&A (special issue for Gaia DR2).        '')')
      write(*,'(''#                                                '')')
      write(*,'(''#                                                '')')
 
      STOP
      end
 
 
 
 
c--------------------------------------------------------
c
c sub_query_gaia.F -
c
c   this is a FORTRAN subroutine that queries the GAIA catalog
c
c--------------------------------------------------------
 
      subroutine sub_query_GAIA_DR3(RA_CEN,DE_CEN,RAD,
     .                              ra_o,  de_o, er_o, ed_o,
     .                              mg_o, emg_o,
     .                              mb_o, emb_o,
     .                              mr_o, emr_o,
     .                            mura_o, emura_o,
     .                            mude_o, emude_o,
     .                              pi_o, epi_o,
     .                              id_o,
     .                             b85_o,
     .                             Os,DIAG)
      implicit none
 
      real*8  RA_CEN
      real*8  DE_CEN
      real*4  RAD
 
      real*8      ra_o(_GMAX_)
      real*8      de_o(_GMAX_)
      real*8      er_o(_GMAX_)
      real*8      ed_o(_GMAX_)
      real*8      mg_o(_GMAX_)
      real*8     emg_o(_GMAX_)
      real*8      mb_o(_GMAX_)
      real*8     emb_o(_GMAX_)
      real*8      mr_o(_GMAX_)
      real*8     emr_o(_GMAX_)
      real*8    mura_o(_GMAX_)
      real*8   emura_o(_GMAX_)
      real*8    mude_o(_GMAX_)
      real*8   emude_o(_GMAX_)
      real*8      pi_o(_GMAX_)
      real*8     epi_o(_GMAX_)
      integer*8   id_o(_GMAX_)
      byte       b85_o(85,_GMAX_)
 
      integer Os
      logical DIAG
 
      integer  n_o(_GMAX_)
      real    dr_o(_GMAX_)
      real    dd_o(_GMAX_)
 
      real*8 rd2x, rd2y
      real dra, dde, ddd
 
c---------------------------------------
c
 
      integer i0
      character*80 stringg
      character*80 STRING
 
      integer nstar
 
      real*8    RA       ! RA in degrees
      real*8    DE       ! De in degrees
      real*8   eRA       ! error in RA (deg)
      real*8   eDE       ! error in DE (deg)
 
      real*8  GMAG       ! GMAG (in mag)
      real*8 eGMAG       ! error in GMAG
      real*8  BMAG       ! BMAG
      real*8 eBMAG       ! error in BMAG
      real*8  RMAG       ! RMAG
      real*8 eRMAG       ! error in RMAG
 
      real*8  muRA       ! pmRA (in arcsec/year)
      real*8 emuRA       ! error in pmRA
      real*8  muDE       ! pmDE
      real*8 emuDE       ! error in pmDE
      real*8    PI       ! parallax
      real*8   ePI       ! error in parralax
 
      integer*8 LAB_I8   ! GAIA label
      integer*2  APARU   ! ASTROM params
 
      integer*8 NVIS     ! number of visibilities
      integer*8 NMAT     !
      real*8    IGOF     !
      real*8    IXSN     !
      integer*8 NBAD     ! number of bad measts
      integer*8 NGUD     ! number of good measts
      integer*8 UFLG     ! use flag
      integer*8 GFLG     ! good flag
      integer*8 RFLG     ! is there a RV?
      real*8      RV     ! rv
      real*8     ERV     ! rv error
 
 
c-----------------------------------
c
c
c     byte buff50_100(50,100)
c     byte buff5000(5000)
c     equivalence(buff5000,buff50_100)
c     byte locrec_b(50)
 
      byte buff85_100(85,100)
      byte buff8500(8500)
      equivalence(buff8500,buff85_100)
 
      byte locrec_b85(85)
 
      integer N
 
c-------------------------------------
c
c this map tells you what the first record
c is for any point in the sky
c
      integer*4 NMAP_IJ(3601,1801)
      byte      NMAP_B(25941604)
      equivalence(NMAP_IJ,NMAP_B)
      common  /NMAP_IJ_/NMAP_IJ
 
      integer   Qs
c     byte      B_12L(12,7224631)
c     integer*8 ID_Q(7224631)
c     real*4    RV_Q(7224631)
c     real*4   eRV_Q(7224631)
c     common  /RV_Q_/ID_Q, RV_Q, eRV_Q, Qs
 
      logical FIRST
      data    FIRST / .true. /
      common /FIRST_/FIRST
 
 
 
c---------------------------------------------------------------
 
      integer IIMINu,IIMAXu
      integer IIMIN, IIMAX, IICEN
      integer JJMIN, JJMAX, JJCEN
 
      integer i, j
 
 
      real*8 rminij, rmaxij
      real*8 dminij, dmaxij
 
      integer Ls
      integer B
 
      real*8 RA_MIN, RA_MAX
      real*8 DE_MIN, DE_MAX
 
      integer ios
      integer ii, jj
 
      integer n1, n2
      integer r1, r2
 
c     integer*8 I8, I8_RV, I8_RVE
c     byte      B1_8(8)
c     equivalence (I8,B1_8)
 
      integer*2 bsho
      logical   unselect
      integer   UWHY
      integer*8 r
      integer*8 q
 
      if (FIRST) then
         if (DIAG) then
             write(*,'(''#FIRST = '',l1)') FIRST
             write(*,'(''#open BASE.GDR2 for map image...'')')
             endif
         Qs = 0
         open(99,file='NONE',
     .           action='read',
     .           status='old',
     .           err=900,recl=25941604,form='UNFORMATTED',
     .           access='DIRECT')
         read(99,rec=1,iostat=ios) NMAP_B
         close(99)
         endif
 
      if (DIAG) then
         print*,' '
         print*,'sub_query_GAIA ---  '
         print*,'    RA_CEN: ',RA_CEN
         print*,'    DE_CEN: ',DE_CEN
         print*,'       RAD: ',RAD
         print*,' '
         endif
 
      RA_MIN = RA_CEN - RAD/cos(DE_CEN*3.14159/180)
      RA_MAX = RA_CEN + RAD/cos(DE_CEN*3.14159/180)
 
      if (abs(DE_CEN).gt.90-RAD) then
          RA_MIN =   0.00
          RA_MAX = 360.00
          endif
 
      if (RA_MIN.lt.0) then
          RA_MIN =   0.0
          RA_MAX = 360.0
          endif
 
      if (RA_MAX.gt.360) then
          RA_MIN =   0.0
          RA_MAX = 360.0
          endif
 
      DE_MIN = DE_CEN - RAD
      DE_MAX = DE_CEN + RAD
 
      if (DIAG) then
         write(*,*)
         write(*,'(5x,a12,3f16.7)') '--> RA_CEN: ',RA_CEN
         write(*,'(5x,a12,3f16.7)') '--> DE_CEN: ',DE_CEN
         write(*,'(5x,a12,3f16.7)') '-->    RAD: ',RAD
         write(*,*)
         write(*,'(5x,a12,3f16.7)') '--> RA_MIN: ',RA_MIN,RA_CEN,RA_MAX
         write(*,'(5x,a12,3f16.7)') '--> DE_MIN: ',DE_MIN,DE_CEN,DE_MAX
         write(*,*)
         endif
 
      if (DE_MIN.lt.-90.0) DE_MIN = -90.0
      if (DE_MAX.gt. 90.0) DE_MAX =  90.0
 
      if (DIAG.and.(.false.)) print*,'OPEN BASE.DATA...'
 
      open(99,file='NONE',
     .        action='read',
     .        status='old',
     .        recl=8500,
     .        form='UNFORMATTED',
     .        access='DIRECT')
 
      iicen = 1 + (RA_MIN   )*10 + 0.5
      jjcen = 1 + (DE_MIN+90)*10 + 0.5
 
      iimin = 1 + RA_MIN*10      + 0.5 - 2
      iimax = 1 + RA_MAX*10      + 0.5 + 2
      jjmin = 1 + (DE_MIN+90)*10 + 0.5 - 2
      jjmax = 1 + (DE_MAX+90)*10 + 0.5 + 2
 
      if (iimin.lt.0001) iimin = 0001
      if (iimax.gt.3601) iimax = 3601
      if (jjmin.lt.0001) jjmin = 0001
      if (jjmax.gt.1801) jjmax = 1801
 
      if (DIAG) then
         write(*,*)
         write(*,*)
         write(*,113) 'iminmax: ',iimin,iicen,iimax,
     .                            (iimin-1)*0.1,
     .                            (iicen-1)*0.1,
     .                            (iimax-1)*0.1
         write(*,113) 'jminmax: ',jjmin,jjcen,jjmax,
     .                            (jjmin-1)*0.1-90,
     .                            (jjcen-1)*0.1-90,
     .                            (jjmax-1)*0.1-90
  113    format(1x,a10,1x,3i4,1x,3f15.7)
         write(*,*)
         do jj = jjmin, jjmax
         write(35,*)
         write(*,*)
         do ii = iimin, iimax+1
            write(*,114) ii,jj,nmap_ij(ii,jj),
     .                         nmap_ij(ii,jj)-nmap_ij(max(ii-1,1),jj)
            if (DIAG.and.(.false.))
     .      write(*,114) ii,jj,nmap_ij(ii,jj),
     .                         nmap_ij(ii,jj)-nmap_ij(max(ii-1,1),jj)
  114       format(20x,i4,1x,i4,1x,2i10)
            enddo
            enddo
         print*,'          '
         endif
 
      if (DIAG) then
          open(38,file='query_GAIA_DR2.diag.01.jrd',status='unknown')
          open(39,file='query_GAIA_DR2.diag.02.chk',status='unknown')
          open(40,file='query_GAIA_DR2.diag.03.uns',status='unknown')
          endif
      Ls = 0
      Os = 0
      do j = jjmin, jjmax
         r1 = 1 + int(abs(NMAP_IJ(iimin  ,j))/100)
         r2 = 1 + int(abs(NMAP_IJ(iimax+1,j))/100)
         if (DIAG) write( *,'(6i12)') j,iimin,iimax+1,r1,r2,r2-r1
         if (DIAG) write(36,'(6i12)') j,iimin,iimax+1,r1,r2,r2-r1
         if (r2.lt.r1) stop 'r2.lt.r1 problem'
         do r = r1, r2
            read(99,rec=r+3052,iostat=ios) buff8500
            if (DIAG) write(37,'(20x,i5,1x,i8,1x,i8)') j, r, Ls
            do q = 1, 100
               if (q+(r-1)*100.le.1811709771) then
                   call b50_rdpmugbv(buff85_100(1,q),
     .                              RA,eRA,  DE,eDE,
     .                              PI, ePI, muRA,emuRA,muDE,emuDE,
     .                              GMAG,eGMAG,BMAG,eBMAG,RMAG,eRMAG,
     .                              LAB_I8, APARU)
                   Ls = Ls + 1
                   dra = rd2x(RA,DE,RA_CEN,DE_CEN)
                   dde = rd2y(RA,DE,RA_CEN,DE_CEN)
                   ddd = sqrt(dra**2+dde**2)
  138              format(1x,i4,1x,i9,1x,f9.2,1x,f9.2,1x,f9.2,1x,i19)
                   if (ddd.lt.RAD) then
                       Os = Os + 1
                       if (Os.le._GMAX_) then
                           ra_o(Os)    =  RA
                           de_o(Os)    =  DE
                           er_o(Os)    = eRA*1000.     ! this function returned it in mas
                           ed_o(Os)    = eDE*1000.     ! this function returned it in mas
                           mg_o(Os)    =  GMAG
                           emg_o(Os)   = eGMAG*1000.   ! this function returned it in mmag
                           mb_o(Os)    =  BMAG
                           emb_o(Os)   = eBMAG*1000.   ! this function returned it in mmag
                           mr_o(Os)    =  RMAG
                           emr_o(Os)   = eRMAG*1000.   ! this function returned it in mmag
                           mura_o(Os)  =  muRA
                           emura_o(Os) = emuRA
                           mude_o(Os)  =  muDE
                           emude_o(Os) = emuDE
                           pi_o(Os)    =  PI
                           epi_o(Os)   = ePI
                           id_o(Os)    = LAB_I8
                           n_o(Os)     = q + (r-1)*100
                           do i = 1, 85
                              b85_o(i,Os) = buff85_100(i,q)
                              enddo
                           if (unselect(RA,DE,eRA,eDE,
     .                                  GMAG,eGMAG,BMAG,eBMAG,
     .                                             RMAG,eRMAG,
     .                                  muRA,emuRA,muDE,emuDE,PI,ePI,
     .                                  NVIS, NMAT, IGOF, IXSN,
     .                                  NBAD, NGUD, UFLG, GFLG, RFLG,
     .                                  RV, eRV,UWHY)) then
                                        Os = Os-1
                                        endif
                           endif
                       endif
                   endif
               enddo
            enddo
         enddo
      if (DIAG) then
          write(*,'(''# TOTAL READ Ls: '',i8)') Ls
          write(*,'(''# TOTAL INCL Os: '',i8)') Os
          endif
      close(38)
      close(99)
 
 
      return
 
 900  continue
      print*,' '
      print*,'subroutine sub_query_GAIA_DR3 ... '
      print*,' '
      print*,'File open error '
      print*,'   ios: ', ios
      print*,' '
      print*,'  file:  NONE'
      print*,' '
      stop
      end
 
 
 
 
 
 
      subroutine sub_query_GAIA_DR3_DUMP(RA_MINU,RA_MAXU,
     .                                   DE_MINU,DE_MAXU)
      implicit none
 
      real*8   RA_MINU,    RA_MAXU
      real*8   DE_MINU,    DE_MAXU
 
      integer*8 r, q
 
      byte buff85_100(85,100)
      byte buff8500(8500)
      equivalence(buff8500,buff85_100)
 
      integer ios
      integer iii
 
      integer*4 NMAP_IJ(3601,1801)
      byte      NMAP_B(25941604)
      equivalence(NMAP_IJ,NMAP_B)
      common  /NMAP_IJ_/NMAP_IJ
 
      integer i1, i2, i
      integer j1, j2, j
 
      integer r1_j(1801), r1
      integer r2_j(1801), r2
 
      logical unselect_b85
 
      do j = 0001, 1801
         r1_j(j) = -1
         r2_j(j) = -1
         enddo
 
      open(99,file='NONE',
     .        action='read',
     .        status='old',
     .        iostat=ios,
     .        err=900,recl=25941604,form='UNFORMATTED',
     .        access='DIRECT')
      read(99,rec=1,iostat=ios) NMAP_B
      close(99)
 
      i1 = max(0001,int(RA_MINU*10)-1)
      i2 = min(3601,int(RA_MAXU*10)+1)
      j1 = max(0001,int((DE_MINU+90)*10)-1)
      j2 = min(1801,int((DE_MAXU+90)*10)+1)
 
c
c     print*,'           '
c     print*,'  RA_MINU: ',RA_MINU, i1
c     print*,'  RA_MAXU: ',RA_MAXU, i2
c     print*,'           '
c     print*,'  DE_MINU: ',DE_MINU, j1
c     print*,'  DE_MAXU: ',DE_MAXU, j2
c     print*,'           '
c
 
      open(99,file='NONE',
     .        action='read',
     .        status='old',
     .        recl=8500,
     .        form='UNFORMATTED',
     .        access='DIRECT')
 
      r1 = nmap_ij(3601,1801)
      r2 = nmap_ij(0001,0001)
      do j = j1, j2
         r1_j(j) = nmap_ij(i1,j)/100
         r2_j(j) = nmap_ij(i2,j)/100+1
         if (i2.eq.3601.and.j2.ne.1801)
     .       r2_j(j) = nmap_ij(0001,j+1)/100+1
         if (r1_j(j).le.r1) r1 = r1_j(j)
         if (r2_j(j).ge.r2) r2 = r2_j(j)
         enddo
 
 
c     print*,' '
c     print*,'  r1: ',r1
c     print*,'  r2: ',r2
c     print*,' '
 
      if (r1.lt.       1) r1 =        1
      if (r2.gt.18117098) r2 = 18117098
 
      iii = 0
c     do r = 1, 18117098
      do r = r1, r2
         do j = 0001, 1801
            if (r.ge.r1_j(j).and.r.le.r2_j(j)) goto 1
            enddo
         goto 2
   1     continue
         if (r1.le.1.and.r/18117.ge.iii) then
             write(0,119) r, 100.*r/18117098, iii
 119         format(1x,i10,1x,f8.3,' % done ',1x,i4)
             iii = iii + 1
             endif
         read(99,rec=r+3052,iostat=ios) buff8500
         do q = 1, 100
            if (.not.unselect_b85(buff85_100(1,q))) then
                if (q+(r-1)*100.le.1811709771) then
                    call showfull_b85(buff85_100(1,q),.false.,6)
                    endif
                endif
            enddo
   2     continue
         enddo
 
      close(99)
      close(98)
      return
 
  900 continue
      print*,'                        '
      print*,'subroutine sub_query_GAIA_DR3_DUMP... '
      print*,'  ios: ',ios
      print*,'                         '
      print*,'  likely problem opening  NONE'
      print*,'                         '
      stop
      end
 
 
 
 
 
 
c-----------------------------------------------
c
c This routine will take an (ra,dec) and a
c tangent-point (ra0,dec0) and will compute
c the y coordinate in the tangent plane.  The
c tangent-plane y-axis is directed along +DEC.
c The units for y will be degrees, just lilke
c those for RA and DEC.
c
      real*8 function rd2y(r,d,r0,d0)
      implicit none
 
      real*8 r, d
      real*8 r0,d0
 
      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 yrad
      real   x, y, z
      real   xx,yy,zz
 
      cosra = cos((r-r0)*3.141592654d0/180.0d0)
      sinra = sin((r-r0)*3.141592654d0/180.0d0)
      cosde = cos(d *3.141592654d0/180.0d0)
      sinde = sin(d *3.141592654d0/180.0d0)
      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)
 
      rrrr = sind0*sinde + cosd0*cosde*cosra
 
      yrad = (cosd0*sinde-sind0*cosde*cosra)/rrrr
      rd2y = yrad*180.0d0/3.141592654d0
 
      x  = cosde*cos(r *3.14159/180)
      y  = cosde*sin(r *3.14159/180)
      z  = sinde
      xx = cosd0*cos(r0*3.14159/180)
      yy = cosd0*sin(r0*3.14159/180)
      zz = sind0
 
      if (x*xx + y*yy + z*zz.lt.0) rd2y = 90
 
      return
      end
 
 
c-----------------------------------------------
c
c This routine will take an (ra,dec) and a
c tangent-point (ra0,dec0) and will compute
c the x coordinate in the tangent plane.  The
c tangent-plane x-axis is directed along -RA.
c The units for x will be degrees, just lilke
c those for RA and DEC.
c
      real*8 function rd2x(r,d,r0,d0)
      implicit none
 
      real*8 r, d
      real*8 r0,d0
 
      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 xrad
      real   x, y, z
      real   xx,yy,zz
 
      cosra = cos((r-r0)*3.141592654d0/180.0d0)
      sinra = sin((r-r0)*3.141592654d0/180.0d0)
 
      cosde = cos(d *3.141592654d0/180.0d0)
      sinde = sin(d *3.141592654d0/180.0d0)
      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)
 
      rrrr = sind0*sinde + cosd0*cosde*cosra
 
 
      xrad = cosde*sinra/rrrr
      rd2x = xrad*180.0d0/3.141592654d0
 
      x  = cosde*cos(r *3.14159/180)
      y  = cosde*sin(r *3.14159/180)
      z  = sinde
      xx = cosd0*cos(r0*3.14159/180)
      yy = cosd0*sin(r0*3.14159/180)
      zz = sind0
 
      if (x*xx + y*yy + z*zz.lt.0) rd2x = 90
 
      return
      end
 
 
 
 
 
 
 
c--------------------------------------------------------------------
c
c
      byte function use(r,d,mj,mh,mk,
     .                  rcen,dcen,jcen,hcen,kcen,
     .                       drmin,drmax,
     .                       djmin,djmax,dhmin,dhmax,dkmin,dkmax,
     .                       dxmax,dymax,
     .                  rmin,rmax,dmin,dmax,
     .                  ebmin,ebmax,elmin,elmax,
     .                  gbmin,gbmax,glmin,glmax,
     .                  jmin,jmax,hmin,hmax,kmin,kmax)
      implicit none
 
      real*8 r,d
      real   mj,mh,mk
 
      real*8 rcen, dcen
      real   jcen, hcen, kcen
      real   drmin, drmax
      real   djmin, djmax, dhmin, dhmax, dkmin, dkmax
      real   dxmax, dymax
      real*8  rmin,  rmax,  dmin,  dmax
      real*8 ebmin, ebmax, elmin, elmax
      real*8 gbmin, gbmax, glmin, glmax
      real   jmin, jmax
      real   hmin, hmax
      real   kmin, kmax
 
      real*8 gl, gb
      real*8 el, eb
 
      real*8 rd2x, rd2y
      real*8 dx, dy, dd
 
      if (rmin.gt.0) then
         use = -1
         if (r.lt.rmin) return
         use = -2
         if (r.gt.rmax) return
         use = -3
         if (d.lt.dmin) return
         use = -4
         if (d.gt.dmax) return
         endif
 
      if (jmin.gt.-9.99) then
         use = -13
         if (mj.lt.jmin) return
         use = -14
         if (mj.gt.jmax) return
         endif
 
      if (hmin.gt.-9.99) then
         use = -15
         if (mh.lt.hmin) return
         use = -16
         if (mh.gt.hmax) return
         endif
 
      if (kmin.gt.-9.99) then
         use = -17
         if (mk.lt.kmin) return
         use = -18
         if (mk.gt.kmax) return
         endif
 
      if (rcen.ge.0.00) then
         dx = rd2x(r,d,rcen,dcen)
         dy = rd2y(r,d,rcen,dcen)
         dd = sqrt(dx**2+dy**2)
         use = -19
         if (dd.lt.drmin) return
         use = -20
         if (dd.gt.drmax) return
         use = -21
         if (mj-jcen.gt.djmax) return
         use = -22
         if (mj-jcen.lt.djmin) return
         use = -23
         if (mh-hcen.gt.dhmax) return
         use = -24
         if (mh-hcen.lt.dhmin) return
         use = -25
         if (mk-kcen.gt.dkmax) return
         use = -26
         if (mk-kcen.lt.dkmin) return
         use = -27
         if (abs(dx).gt.dxmax) return
         use = -28
         if (abs(dy).gt.dymax) return
         endif
 
      use = 1
 
      return
      end
 
 
C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.
 
 
      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      REAL A(N)
      REAL X,XX,Z,ZZ,Y
C
      NDEEP=0
      U1=N
      L1=1
      DO I=1,N
    1    ORD(I)=I
         enddo
    2 IF (U1.LE.L1) RETURN
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END
 
 
 
 
c------------------------------------------------
c
c same as rclip, but for double precision numbers
c
      real*8 function dclip(dval,dlo,dhi)
      implicit none
 
      real*8 dval
      real*8 dlo
      real*8 dhi
 
      dclip = dval
      if (dclip.gt.dhi) dclip = dhi
      if (dclip.lt.dlo) dclip = dlo
      if (.not.(dclip.lt.dhi).and.
     .    .not.(dclip.gt.dlo)) dclip = dhi
 
      return
      end
 
 
 
 
c------------------------------------------------
c
c turn a byte into something unsigned we can show
c
      integer*2 function bsho(b)
      implicit none
 
      byte b
      bsho = b
      if (bsho.lt.0) bsho = bsho + 256
 
      return
      end
 
 
 
 
c--------------------------------------------
c
c this function will return the record where
c the star ID is located
c
c
c      integer*4 function r_id(id)
c      implicit none
c
c      integer*8 id
c
c      integer*8 nst1, id1, r1
c      integer*8 nst2, id2, r2
c      integer*8 nstb, idb, rb
c
c      integer*8 NSTs
c
c      integer NIT
c
c      NSTs = 1 692 919 135
c
cc     write(*,116) 0,id
c      nst1 = 1
c      nst2 = NSTs
c
c
c      call nst2idNr(nst1,id1,r1)
c      if (id.eq.id1) then
c          r_id = r1
c          return
c          endif
c
c      call nst2idNr(nst2,id2,r2)
c      if (id.eq.id2) then
c          r_id = r2
c          return
c          endif
c
c      NIT = 0
c    1 continue
c      NIT = NIT + 1
c      if (NIT.gt.50) then
c          print*,' '
c          print*,'Problem with r_id ...   '
c          print*,'Cannot find record for star ID: ',id
c          print*,'After more than 50 binary search moves.'
c          print*,' '
c          stop
c          endif
c      nstb = (nst1+nst2)/2
c      call nst2idNr(nstb,idb,rb)
cc     write(*,116) NIT, id,idb,nstb,rb,
cc    .                     id1,nst1,r1,
cc    .                     id2,nst2,r2
c  116 format(i4,1x,i19,3(10x,i19,1x,i10,1x,i8))
c      if (id.gt.idb) then
c          nst1 = nstb
c          id1  = idb
c          r1   = rb
c          goto 1
c          endif
c      if (id.lt.idb) then
c          nst2 = nstb
c          id2  = idb
c          r2   = rb
c          goto 1
c          endif
cc     write(*,116) NIT, id,idb,nstb,rb
c
c      r_id = rb
c
c      return
c      end
c
c
 
 
 
 
 
      logical function unselect_b85(b_85)
 
      implicit none
 
      byte b_85(85)
 
c--------------------------------------------
 
      real*8    RA       ! RA in degrees
      real*8    DE       ! De in degrees
      real*8   eRA       ! error in RA (deg)
      real*8   eDE       ! error in DE (deg)
 
      real*8  GMAG       ! GMAG (in mag)
      real*8 eGMAG       ! error in GMAG
      real*8  BMAG       ! BMAG
      real*8 eBMAG       ! error in BMAG
      real*8  RMAG       ! RMAG
      real*8 eRMAG       ! error in RMAG
 
      real*8  muRA       ! pmRA (in mas/year)
      real*8 emuRA       ! error in pmRA
      real*8  muDE       ! pmDE (in mas)
      real*8 emuDE       ! error in pmDE
      real*8    PI       ! parallax (in mas)
      real*8   ePI       ! error in parralax
 
      integer*8 SID
      integer*2 apars
 
      integer   UWHY     ! why to unselect
 
      real*8   RA_MINU,    RA_MAXU
      real*8   DE_MINU,    DE_MAXU
      real*8  ERA_MINU,   ERA_MAXU
      real*8  EDE_MINU,   EDE_MAXU
      real*8   MG_MINU,    MG_MAXU
      real*8  EMG_MINU,   EMG_MAXU
      real*8   MB_MINU,    MB_MAXU
      real*8  EMB_MINU,   EMB_MAXU
      real*8   MR_MINU,    MR_MAXU
      real*8   EMR_MINU,  EMR_MAXU
      real*8  MURA_MINU, MURA_MAXU
      real*8 EMURA_MINU,EMURA_MAXU
      real*8  MUDE_MINU, MUDE_MAXU
      real*8 EMUDE_MINU,EMUDE_MAXU
      real*8    PI_MINU,   PI_MAXU
      real*8   EPI_MINU,  EPI_MAXU
      real*8  NVIS_MINU, NVIS_MAXU
      real*8  NMAT_MINU, NMAT_MAXU
      real*8  IGOF_MINU, IGOF_MAXU
      real*8  IXSN_MINU, IXSN_MAXU
      real*8  NBAD_MINU, NBAD_MAXU
      real*8  NGUD_MINU, NGUD_MAXU
      real*8  UFLG_MINU, UFLG_MAXU
      real*8  GFLG_MINU, GFLG_MAXU
      real*8  RFLG_MINU, RFLG_MAXU
      real*8    RV_MINU,   RV_MAXU
      real*8   ERV_MINU,  ERV_MAXU
      real*8   BmR_MINU,  BmR_MAXU
      real*8   BmG_MINU,  BmG_MAXU
      real*8   GmR_MINU,  GmR_MAXU
 
      real*8 BmG, BmR, GmR
 
      common / LIMITS_ /  RA_MINU,    RA_MAXU,
     .                    DE_MINU,    DE_MAXU,
     .                   ERA_MINU,   ERA_MAXU,
     .                   EDE_MINU,   EDE_MAXU,
     .                    MG_MINU,    MG_MAXU,
     .                   EMG_MINU,   EMG_MAXU,
     .                    MB_MINU,    MB_MAXU,
     .                   EMB_MINU,   EMB_MAXU,
     .                    MR_MINU,    MR_MAXU,
     .                    EMR_MINU,  EMR_MAXU,
     .                   MURA_MINU, MURA_MAXU,
     .                  EMURA_MINU,EMURA_MAXU,
     .                   MUDE_MINU, MUDE_MAXU,
     .                  EMUDE_MINU,EMUDE_MAXU,
     .                     PI_MINU,   PI_MAXU,
     .                    EPI_MINU,  EPI_MAXU,
     .                   NVIS_MINU, NVIS_MAXU,
     .                   NMAT_MINU, NMAT_MAXU,
     .                   IGOF_MINU, IGOF_MAXU,
     .                   IXSN_MINU, IXSN_MAXU,
     .                   NBAD_MINU, NBAD_MAXU,
     .                   NGUD_MINU, NGUD_MAXU,
     .                   UFLG_MINU, UFLG_MAXU,
     .                   GFLG_MINU, GFLG_MAXU,
     .                   RFLG_MINU, RFLG_MAXU,
     .                     RV_MINU,   RV_MAXU,
     .                    ERV_MINU,  ERV_MAXU,
     .                    BmR_MINU,  BmR_MAXU,
     .                    BmG_MINU,  BmG_MAXU,
     .                    GmR_MINU,  GmR_MAXU
 
      unselect_b85 = .false.
 
      call b50_rdpmugbv(b_85,
     .                  RA,eRA,
     .                  DE,eDE,
     .                  PI,ePI,
     .                  muRA,emuRA,
     .                  muDE,emuDE,
     .                  GMAG,eGMAG,
     .                  BMAG,eBMAG,
     .                  RMAG,erMAG,
     .                  SID,
     .                  apars)
 
      BmG = BMAG - GMAG
      BmR = BMAG - RMAG
      GmR = GMAG - RMAG
 
      UWHY = 0
      if (   RA.lt.   RA_MINU) UWHY =   1
      if (   RA.gt.   RA_MAXU) UWHY =   2
      if (   DE.lt.   DE_MINU) UWHY =   3
      if (   DE.gt.   DE_MAXU) UWHY =   4
      if (  eRA.lt.  ERA_MINU) UWHY =   5
      if (  eRA.gt.  ERA_MAXU) UWHY =   6
      if (  eDE.lt.  EDE_MINU) UWHY =   7
      if (  eDE.gt.  EDE_MAXU) UWHY =   8
      if ( GMAG.lt.   MG_MINU) UWHY =   9
      if ( GMAG.gt.   MG_MAXU) UWHY =  10
      if (eGMAG.lt.  EMG_MINU) UWHY =  11
      if (eGMAG.gt.  EMG_MAXU) UWHY =  12
      if ( BMAG.lt.   MB_MINU) UWHY =  13
      if ( BMAG.gt.   MB_MAXU) UWHY =  14
      if (eBMAG.lt.  EMB_MINU) UWHY =  15
      if (eBMAG.gt.  EMB_MAXU) UWHY =  16
      if ( RMAG.lt.   MR_MINU) UWHY =  17
      if ( RMAG.gt.   MR_MAXU) UWHY =  18
      if (eRMAG.lt.  EMR_MINU) UWHY =  19
      if (eRMAG.gt.  EMR_MAXU) UWHY =  20
      if ( muRA.lt. MURA_MINU) UWHY =  21
      if ( muRA.gt. MURA_MAXU) UWHY =  22
      if (emuRA.lt.EMURA_MINU) UWHY =  23
      if (emuRA.gt.EMURA_MAXU) UWHY =  24
      if ( muDE.lt. MUDE_MINU) UWHY =  25
      if ( muDE.gt. MUDE_MAXU) UWHY =  26
      if (emuDE.lt.EMUDE_MINU) UWHY =  27
      if (emuDE.gt.EMUDE_MAXU) UWHY =  28
      if (   PI.lt.   PI_MINU) UWHY =  29
      if (   PI.gt.   PI_MAXU) UWHY =  30
      if (  ePI.lt.  EPI_MINU) UWHY =  31
      if (  ePI.gt.  EPI_MAXU) UWHY =  32
      if (  BmG.lt.  BmG_MINU) UWHY =  55
      if (  BmG.gt.  BmG_MAXU) UWHY =  56
      if (  BmR.lt.  BmR_MINU) UWHY =  57
      if (  BmR.gt.  BmR_MAXU) UWHY =  58
      if (  GmR.lt.  GmR_MINU) UWHY =  59
      if (  GmR.gt.  GmR_MAXU) UWHY =  60
 
      unselect_b85 = UWHY.ne.0
 
      return
      end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
      logical function unselect(RA,DE,eRA,eDE,
     .                          GMAG,eGMAG,BMAG,eBMAG,RMAG,eRMAG,
     .                          muRA,emuRA,muDE,emuDE,PI,ePI,
     .                          NVIS, NMAT, IGOF, IXSN,
     .                          NBAD, NGUD, UFLG, GFLG, RFLG,
     .                          RV, ERV, UWHY)
      implicit none
 
      real*8    RA       ! RA in degrees
      real*8    DE       ! De in degrees
      real*8   eRA       ! error in RA (deg)
      real*8   eDE       ! error in DE (deg)
 
      real*8  GMAG       ! GMAG (in mag)
      real*8 eGMAG       ! error in GMAG
      real*8  BMAG       ! BMAG
      real*8 eBMAG       ! error in BMAG
      real*8  RMAG       ! RMAG
      real*8 eRMAG       ! error in RMAG
 
      real*8  muRA       ! pmRA (in mas/year)
      real*8 emuRA       ! error in pmRA
      real*8  muDE       ! pmDE (in mas)
      real*8 emuDE       ! error in pmDE
      real*8    PI       ! parallax (in mas)
      real*8   ePI       ! error in parralax
 
      integer*8 NVIS     ! number of visibilities
      integer*8 NMAT     !
      real*8    IGOF     !
      real*8    IXSN     !
      integer*8 NBAD     ! number of bad measts
      integer*8 NGUD     ! number of good measts
      integer*8 UFLG     ! use flag
      integer*8 GFLG     ! good flag
      integer*8 RFLG     ! is there a RV?
      real*8      RV     ! rv
      real*8     ERV     ! rv error
 
      integer   UWHY     ! why to unselect
 
      real*8   RA_MINU,    RA_MAXU
      real*8   DE_MINU,    DE_MAXU
      real*8  ERA_MINU,   ERA_MAXU
      real*8  EDE_MINU,   EDE_MAXU
      real*8   MG_MINU,    MG_MAXU
      real*8  EMG_MINU,   EMG_MAXU
      real*8   MB_MINU,    MB_MAXU
      real*8  EMB_MINU,   EMB_MAXU
      real*8   MR_MINU,    MR_MAXU
      real*8   EMR_MINU,  EMR_MAXU
      real*8  MURA_MINU, MURA_MAXU
      real*8 EMURA_MINU,EMURA_MAXU
      real*8  MUDE_MINU, MUDE_MAXU
      real*8 EMUDE_MINU,EMUDE_MAXU
      real*8    PI_MINU,   PI_MAXU
      real*8   EPI_MINU,  EPI_MAXU
      real*8  NVIS_MINU, NVIS_MAXU
      real*8  NMAT_MINU, NMAT_MAXU
      real*8  IGOF_MINU, IGOF_MAXU
      real*8  IXSN_MINU, IXSN_MAXU
      real*8  NBAD_MINU, NBAD_MAXU
      real*8  NGUD_MINU, NGUD_MAXU
      real*8  UFLG_MINU, UFLG_MAXU
      real*8  GFLG_MINU, GFLG_MAXU
      real*8  RFLG_MINU, RFLG_MAXU
      real*8    RV_MINU,   RV_MAXU
      real*8   ERV_MINU,  ERV_MAXU
      real*8   BmR_MINU,  BmR_MAXU
      real*8   BmG_MINU,  BmG_MAXU
      real*8   GmR_MINU,  GmR_MAXU
 
      real*8 BmG, BmR, GmR
 
      common / LIMITS_ /  RA_MINU,    RA_MAXU,
     .                    DE_MINU,    DE_MAXU,
     .                   ERA_MINU,   ERA_MAXU,
     .                   EDE_MINU,   EDE_MAXU,
     .                    MG_MINU,    MG_MAXU,
     .                   EMG_MINU,   EMG_MAXU,
     .                    MB_MINU,    MB_MAXU,
     .                   EMB_MINU,   EMB_MAXU,
     .                    MR_MINU,    MR_MAXU,
     .                    EMR_MINU,  EMR_MAXU,
     .                   MURA_MINU, MURA_MAXU,
     .                  EMURA_MINU,EMURA_MAXU,
     .                   MUDE_MINU, MUDE_MAXU,
     .                  EMUDE_MINU,EMUDE_MAXU,
     .                     PI_MINU,   PI_MAXU,
     .                    EPI_MINU,  EPI_MAXU,
     .                   NVIS_MINU, NVIS_MAXU,
     .                   NMAT_MINU, NMAT_MAXU,
     .                   IGOF_MINU, IGOF_MAXU,
     .                   IXSN_MINU, IXSN_MAXU,
     .                   NBAD_MINU, NBAD_MAXU,
     .                   NGUD_MINU, NGUD_MAXU,
     .                   UFLG_MINU, UFLG_MAXU,
     .                   GFLG_MINU, GFLG_MAXU,
     .                   RFLG_MINU, RFLG_MAXU,
     .                     RV_MINU,   RV_MAXU,
     .                    ERV_MINU,  ERV_MAXU,
     .                    BmR_MINU,  BmR_MAXU,
     .                    BmG_MINU,  BmG_MAXU,
     .                    GmR_MINU,  GmR_MAXU
 
      unselect = .false.
 
      BmG = BMAG - GMAG
      BmR = BMAG - RMAG
      GmR = GMAG - RMAG
 
      UWHY = 0
      if (   RA.lt.   RA_MINU) UWHY =   1
      if (   RA.gt.   RA_MAXU) UWHY =   2
      if (   DE.lt.   DE_MINU) UWHY =   3
      if (   DE.gt.   DE_MAXU) UWHY =   4
      if (  eRA.lt.  ERA_MINU) UWHY =   5
      if (  eRA.gt.  ERA_MAXU) UWHY =   6
      if (  eDE.lt.  EDE_MINU) UWHY =   7
      if (  eDE.gt.  EDE_MAXU) UWHY =   8
      if ( GMAG.lt.   MG_MINU) UWHY =   9
      if ( GMAG.gt.   MG_MAXU) UWHY =  10
      if (eGMAG.lt.  EMG_MINU) UWHY =  11
      if (eGMAG.gt.  EMG_MAXU) UWHY =  12
      if ( BMAG.lt.   MB_MINU) UWHY =  13
      if ( BMAG.gt.   MB_MAXU) UWHY =  14
      if (eBMAG.lt.  EMB_MINU) UWHY =  15
      if (eBMAG.gt.  EMB_MAXU) UWHY =  16
      if ( RMAG.lt.   MR_MINU) UWHY =  17
      if ( RMAG.gt.   MR_MAXU) UWHY =  18
      if (eRMAG.lt.  EMR_MINU) UWHY =  19
      if (eRMAG.gt.  EMR_MAXU) UWHY =  20
      if ( muRA.lt. MURA_MINU) UWHY =  21
      if ( muRA.gt. MURA_MAXU) UWHY =  22
      if (emuRA.lt.EMURA_MINU) UWHY =  23
      if (emuRA.gt.EMURA_MAXU) UWHY =  24
      if ( muDE.lt. MUDE_MINU) UWHY =  25
      if ( muDE.gt. MUDE_MAXU) UWHY =  26
      if (emuDE.lt.EMUDE_MINU) UWHY =  27
      if (emuDE.gt.EMUDE_MAXU) UWHY =  28
      if (   PI.lt.   PI_MINU) UWHY =  29
      if (   PI.gt.   PI_MAXU) UWHY =  30
      if (  ePI.lt.  EPI_MINU) UWHY =  31
      if (  ePI.gt.  EPI_MAXU) UWHY =  32
c     if ( NVIS.lt. NVIS_MINU) UWHY =  33
c     if ( NVIS.gt. NVIS_MAXU) UWHY =  34
c     if ( NMAT.lt. NMAT_MINU) UWHY =  35
c     if ( NMAT.gt. NMAT_MAXU) UWHY =  36
c     if ( IGOF.lt. IGOF_MINU) UWHY =  37
c     if ( IGOF.gt. IGOF_MAXU) UWHY =  38
c     if ( IXSN.lt. IXSN_MINU) UWHY =  39
c     if ( IXSN.gt. IXSN_MAXU) UWHY =  40
c     if ( NBAD.lt. NBAD_MINU) UWHY =  41
c     if ( NBAD.gt. NBAD_MAXU) UWHY =  42
c     if ( NGUD.lt. NGUD_MINU) UWHY =  43
c     if ( NGUD.gt. NGUD_MAXU) UWHY =  44
c     if ( UFLG.lt. UFLG_MINU) UWHY =  45
c     if ( UFLG.gt. UFLG_MAXU) UWHY =  46
c     if ( GFLG.lt. GFLG_MINU) UWHY =  47
c     if ( GFLG.gt. GFLG_MAXU) UWHY =  48
c     if ( RFLG.lt. RFLG_MINU) UWHY =  49
c     if ( RFLG.gt. RFLG_MAXU) UWHY =  50
c     if (   RV.lt.   RV_MINU) UWHY =  51
c     if (   RV.gt.   RV_MAXU) UWHY =  52
c     if (  eRV.lt.  ERV_MINU) UWHY =  53
c     if (  eRV.gt.  ERV_MAXU) UWHY =  54
      if (  BmG.lt.  BmG_MINU) UWHY =  55
      if (  BmG.gt.  BmG_MAXU) UWHY =  56
      if (  BmR.lt.  BmR_MINU) UWHY =  57
      if (  BmR.gt.  BmR_MAXU) UWHY =  58
      if (  GmR.lt.  GmR_MINU) UWHY =  59
      if (  GmR.gt.  GmR_MAXU) UWHY =  60
 
      unselect = UWHY.ne.0
 
      return
      end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
      subroutine show_delimeters
      implicit none
 
      real*8   RA_MINU,    RA_MAXU
      real*8   DE_MINU,    DE_MAXU
      real*8  ERA_MINU,   ERA_MAXU
      real*8  EDE_MINU,   EDE_MAXU
      real*8   MG_MINU,    MG_MAXU
      real*8  EMG_MINU,   EMG_MAXU
      real*8   MB_MINU,    MB_MAXU
      real*8  EMB_MINU,   EMB_MAXU
      real*8   MR_MINU,    MR_MAXU
      real*8   EMR_MINU,  EMR_MAXU
      real*8  MURA_MINU, MURA_MAXU
      real*8 EMURA_MINU,EMURA_MAXU
      real*8  MUDE_MINU, MUDE_MAXU
      real*8 EMUDE_MINU,EMUDE_MAXU
      real*8    PI_MINU,   PI_MAXU
      real*8   EPI_MINU,  EPI_MAXU
      real*8  NVIS_MINU, NVIS_MAXU
      real*8  NMAT_MINU, NMAT_MAXU
      real*8  IGOF_MINU, IGOF_MAXU
      real*8  IXSN_MINU, IXSN_MAXU
      real*8  NBAD_MINU, NBAD_MAXU
      real*8  NGUD_MINU, NGUD_MAXU
      real*8  UFLG_MINU, UFLG_MAXU
      real*8  GFLG_MINU, GFLG_MAXU
      real*8  RFLG_MINU, RFLG_MAXU
      real*8    RV_MINU,   RV_MAXU
      real*8   ERV_MINU,  ERV_MAXU
      real*8   BmR_MINU,  BmR_MAXU
      real*8   BmG_MINU,  BmG_MAXU
      real*8   GmR_MINU,  GmR_MAXU
 
 
      common / LIMITS_ /  RA_MINU,    RA_MAXU,
     .                    DE_MINU,    DE_MAXU,
     .                   ERA_MINU,   ERA_MAXU,
     .                   EDE_MINU,   EDE_MAXU,
     .                    MG_MINU,    MG_MAXU,
     .                   EMG_MINU,   EMG_MAXU,
     .                    MB_MINU,    MB_MAXU,
     .                   EMB_MINU,   EMB_MAXU,
     .                    MR_MINU,    MR_MAXU,
     .                    EMR_MINU,  EMR_MAXU,
     .                   MURA_MINU, MURA_MAXU,
     .                  EMURA_MINU,EMURA_MAXU,
     .                   MUDE_MINU, MUDE_MAXU,
     .                  EMUDE_MINU,EMUDE_MAXU,
     .                     PI_MINU,   PI_MAXU,
     .                    EPI_MINU,  EPI_MAXU,
     .                   NVIS_MINU, NVIS_MAXU,
     .                   NMAT_MINU, NMAT_MAXU,
     .                   IGOF_MINU, IGOF_MAXU,
     .                   IXSN_MINU, IXSN_MAXU,
     .                   NBAD_MINU, NBAD_MAXU,
     .                   NGUD_MINU, NGUD_MAXU,
     .                   UFLG_MINU, UFLG_MAXU,
     .                   GFLG_MINU, GFLG_MAXU,
     .                   RFLG_MINU, RFLG_MAXU,
     .                     RV_MINU,   RV_MAXU,
     .                    ERV_MINU,  ERV_MAXU,
     .                    BmR_MINU,  BmR_MAXU,
     .                    BmG_MINU,  BmG_MAXU,
     .                    GmR_MINU,  GmR_MAXU
 
      write(*,'(''#'')')
      if (   RA_MINU.gt.-9999.or.   RA_MAXU.lt.9999)
     .                       write(*,337) '   RA_MINU',   RA_MINU,
     .                                    '   RA_MAXU',   RA_MAXU
      if (   DE_MINU.gt.-9999.or.   DE_MAXU.lt.9999)
     .                       write(*,337) '   DE_MINU',   DE_MINU,
     .                                    '   DE_MAXU',   DE_MAXU
      if (  ERA_MINU.gt.-9999.or.  ERA_MAXU.lt.9999)
     .                       write(*,337) '  ERA_MINU',  ERA_MINU,
     .                                    '  ERA_MAXU',  ERA_MAXU
      if (  EDE_MINU.gt.-9999.or.  EDE_MAXU.lt.9999)
     .                       write(*,337) '  EDE_MINU',  EDE_MINU,
     .                                    '  EDE_MAXU',  EDE_MAXU
      if (   MG_MINU.gt.-9999.or.   MG_MAXU.lt.9999)
     .                       write(*,337) '   MG_MINU',   MG_MINU,
     .                                    '   MG_MAXU',   MG_MAXU
      if (  EMG_MINU.gt.-9999.or.  EMG_MAXU.lt.9999)
     .                       write(*,337) '  EMG_MINU',  EMG_MINU,
     .                                    '  EMG_MAXU',  EMG_MAXU
      if (   MB_MINU.gt.-9999.or.   MB_MAXU.lt.9999)
     .                       write(*,337) '   MB_MINU',   MB_MINU,
     .                                    '   MB_MAXU',   MB_MAXU
      if (  EMB_MINU.gt.-9999.or.  EMB_MAXU.lt.9999)
     .                       write(*,337) '  EMB_MINU',  EMB_MINU,
     .                                    '  EMB_MAXU',  EMB_MAXU
      if (   MR_MINU.gt.-9999.or.   MR_MAXU.lt.9999)
     .                       write(*,337) '   MR_MINU',   MR_MINU,
     .                                    '   MR_MAXU',   MR_MAXU
      if (  EMR_MINU.gt.-9999.or.  EMR_MAXU.lt.9999)
     .                       write(*,337) '  EMR_MINU',  EMR_MINU,
     .                                    '  EMR_MAXU',  EMR_MAXU
      if ( MURA_MINU.gt.-9999.or. MURA_MAXU.lt.9999)
     .                       write(*,337) ' MURA_MINU', MURA_MINU,
     .                                    ' MURA_MAXU', MURA_MAXU
      if (EMURA_MINU.gt.-9999.or.EMURA_MAXU.lt.9999)
     .                       write(*,337) 'EMURA_MINU',EMURA_MINU,
     .                                    'EMURA_MAXU',EMURA_MAXU
      if ( MUDE_MINU.gt.-9999.or. MUDE_MAXU.lt.9999)
     .                       write(*,337) ' MUDE_MINU', MUDE_MINU,
     .                                    ' MUDE_MAXU', MUDE_MAXU
      if (EMUDE_MINU.gt.-9999.or.EMUDE_MAXU.lt.9999)
     .                       write(*,337) 'EMUDE_MINU',EMUDE_MINU,
     .                                    'EMUDE_MAXU',EMUDE_MAXU
      if (   PI_MINU.gt.-9999.or.   PI_MAXU.lt.9999)
     .                       write(*,337) '   PI_MINU',   PI_MINU,
     .                                    '   PI_MAXU',   PI_MAXU
      if (  EPI_MINU.gt.-9999.or.  EPI_MAXU.lt.9999)
     .                       write(*,337) '  EPI_MINU',  EPI_MINU,
     .                                    '  EPI_MAXU',  EPI_MAXU
      if ( NVIS_MINU.gt.-9999.or. NVIS_MAXU.lt.9999)
     .                       write(*,337) ' NVIS_MINU', NVIS_MINU,
     .                                    ' NVIS_MAXU', NVIS_MAXU
      if ( NMAT_MINU.gt.-9999.or. NMAT_MAXU.lt.9999)
     .                       write(*,337) ' NMAT_MINU', NMAT_MINU,
     .                                    ' NMAT_MAXU', NMAT_MAXU
      if ( IGOF_MINU.gt.-9999.or. IGOF_MAXU.lt.9999)
     .                       write(*,337) ' IGOF_MINU', IGOF_MINU,
     .                                    ' IGOF_MAXU', IGOF_MAXU
      if ( IXSN_MINU.gt.-9999.or. IXSN_MAXU.lt.9999)
     .                       write(*,337) ' IXSN_MINU', IXSN_MINU,
     .                                    ' IXSN_MAXU', IXSN_MAXU
      if ( NBAD_MINU.gt.-9999.or. NBAD_MAXU.lt.9999)
     .                       write(*,337) ' NBAD_MINU', NBAD_MINU,
     .                                    ' NBAD_MAXU', NBAD_MAXU
      if ( NGUD_MINU.gt.-9999.or. NGUD_MAXU.lt.9999)
     .                       write(*,337) ' NGUD_MINU', NGUD_MINU,
     .                                    ' NGUD_MAXU', NGUD_MAXU
      if ( UFLG_MINU.gt.-9999.or. UFLG_MAXU.lt.9999)
     .                       write(*,337) ' UFLG_MINU', UFLG_MINU,
     .                                    ' UFLG_MAXU', UFLG_MAXU
      if ( GFLG_MINU.gt.-9999.or. GFLG_MAXU.lt.9999)
     .                       write(*,337) ' GFLG_MINU', GFLG_MINU,
     .                                    ' GFLG_MAXU', GFLG_MAXU
      if ( RFLG_MINU.gt.-9999.or. RFLG_MAXU.lt.9999)
     .                       write(*,337) ' RFLG_MINU', RFLG_MINU,
     .                                    ' RFLG_MAXU', RFLG_MAXU
      if (   RV_MINU.gt.-9999.or.   RV_MAXU.lt.9999)
     .                       write(*,337) '   RV_MINU',   RV_MINU,
     .                                    '   RV_MAXU',   RV_MAXU
      if (  ERV_MINU.gt.-9999.or.  ERV_MAXU.lt.9999)
     .                       write(*,337) '  ERV_MINU',  ERV_MINU,
     .                                    '  ERV_MAXU',  ERV_MAXU
      if (  BmR_MINU.gt.-9999.or.  BmR_MAXU.lt.9999)
     .                       write(*,337) '  BmR_MINU',  BmR_MINU,
     .                                    '  BmR_MAXU',  BmR_MAXU
      if (  BmG_MINU.gt.-9999.or.  BmG_MAXU.lt.9999)
     .                       write(*,337) '  BmG_MINU',  BmG_MINU,
     .                                    '  BmG_MAXU',  BmG_MAXU
      if (  GmR_MINU.gt.-9999.or.  GmR_MAXU.lt.9999)
     .                       write(*,337) '  GmR_MINU',  GmR_MINU,
     .                                    '  GmR_MAXU',  GmR_MAXU
 337  format('#  RESTRICTION: ',a10,1x,f10.5,1x,a10,f10.5)
      write(*,'(''#'')')
 
      return
      end
 
 
 
 
c-----------------------------------
c
c from ../WCS/sub_calcpi.f
c
c ---> gifted by David Kaplan
c ---> probably isn't the most accurate, but it is simple
c
cfunction [pi_ra,pi_dec]=calcpi(date,ra_rad,dec_rad)
c% function [pi_ra,pi_dec]=calcpi(date,ra_rad,dec_rad)
c
c
c% MATLAB routine to calculate the offsets due to parallax, at a distance of 1 kpc
c% date is in fractional years
c% ra_rad, dec_rad are RA, Dec in radians
c% returns parallactic shift in ra, dec in mas (e.g. 1 kpc)
c% zero-point offset is arbitrary
c
c  time=(date-2000)*365.25;
c  L=(280.466+0.9856474*time);
c  g_rad=(357.528+0.9856003*time)*pi/180;
c  lambda=(L+1.915*sin(g_rad)+.02*sin(2*g_rad));
c  ll_rad=lambda*pi/180;
c  e_rad=(23.44-0.0000004*time)*pi/180;
c  R_au=1.00014 - 0.01671*cos(g_rad) - 0.00014*cos(2.0*g_rad);
c  X = -1 * R_au .* cos(ll_rad);
c  Y = -1 * R_au .* cos(e_rad) .* sin (ll_rad);
c  Z = -1 * R_au .* sin(e_rad) .* sin (ll_rad);
c
c  pi_ra  = X .* sin(ra_rad) - Y.*cos(ra_rad);
c  pi_dec = X .* cos(ra_rad) .* sin(dec_rad)+ Y .* sin(ra_rad) .* sin(dec_rad);
c  pi_dec = pi_dec - Z .* cos(dec_rad);
c
      subroutine sub_calcpi(ra_deg,de_deg,date,par_ra,par_de)
      implicit none
 
      real*8 ra_deg
      real*8 de_deg
      real*8 date
 
      real*8 par_ra
      real*8 par_de
 
      integer i
 
      real*8 time
      real*8 l
      real*8 g_rad
      real*8 lambda
      real*8 ll_rad
      real*8  e_rad
      real*8  R_au
      real*8 x, y, z
      real*8 pi_ra
      real*8 pi_de
 
      real*8 ra_rad
      real*8 de_rad
 
      ra_rad =  ra_deg/180*3.1415927
      de_rad =  de_deg/180*3.1415927
 
      time   = (date-2000)*365.25
      L      = (280.466+0.9856474*time)
      g_rad  = (357.528+0.9856003*time)*3.14159/180
      lambda = (L+1.915*sin(g_rad)+.02*sin(2*g_rad))
      ll_rad = lambda*3.14159/180
      e_rad  = (23.44-0.0000004*time)*3.14159/180
      R_au   = 1.00014 - 0.01671*cos(g_rad) - 0.00014*cos(2.0*g_rad)
      X = -1 * R_au * cos(ll_rad)
      Y = -1 * R_au * cos(e_rad) * sin (ll_rad)
      Z = -1 * R_au * sin(e_rad) * sin (ll_rad)
 
      pi_ra = X * sin(ra_rad) - Y*cos(ra_rad)
      pi_de = X * cos(ra_rad)*sin(de_rad) +
     .        Y * sin(ra_rad)*sin(de_rad)
      pi_de = pi_de - Z * cos(de_rad)
 
      par_ra = pi_ra
      par_de = pi_de
 
      return
      end
 
 
 
 
c-----------------------------------------------------------
c
c from USNO
c
c rad_i, ded_i   --- RA and Dec in degrees
c JDi            --- the input julian date
c
c par_ra, par_de --- the output parallax, in mas
c
c
      subroutine find_pifac_USNO(rad_i,ded_i,rdat_i,parr_o,pard_o)
      implicit none
 
      real*8  rad_i ! deg
      real*8  ded_i ! deg
      real*8 rdat_i ! deg
      real*8 parr_o
      real*8 pard_o
 
      real*8  JDi
      integer E
      real*8  rai
      real*8  dei
      real*8  jd
      real*8  pmra
      real*8  pmde
      real*8  parlax
      real*8  radvel
 
      real*8 tjd
      real*8 rao_pi0, deo_pi0
      real*8 rao_piU, deo_piU
 
      real*8 rd2x, rd2y
 
      integer i
 
 
c     rai =   5.65833378/360*24
c     dei = -72.06497192
 
      rai = rad_i/360*24
      dei = ded_i
 
      JDi = 2451543.5 + (rdat_i-2000.000d0)*365.25
 
      pmra = 0.0d0
      pmde = 0.0d0
 
      parlax = 1.0d0
      radvel = 0.0d0
 
      TJD = JDi
      CALL ASSTAR (TJD,0,RAI,DEI,PMRA,PMDE,0.0d0 ,RADVEL, RAo_pi0,
     .                                                    DEo_pi0)
      CALL ASSTAR (TJD,0,RAI,DEI,PMRA,PMDE,PARLAX,RADVEL, RAo_piU,
     .                                                    DEo_piU)
      parr_o =  rd2x(RAo_piU*15.0d0,DEo_piU,
     .               RAo_pi0*15.0d0,DEo_pi0)*60*60*1000
      pard_o =  rd2y(RAo_piU*15.0d0,DEo_piU,
     .               RAo_pi0*15.0d0,DEo_pi0)*60*60*1000
 
      return
      end
 
 
 
 
c*********************************************
c****
c**** #include "NOVAS_F3.1.f"
c****
c*********************************************
 
*  NOVAS FORTRAN VERS F3.1 of 2011 MARCH 21
 
************************************************************************
*                                                                      *
*                              N O V A S                               *
*           NAVAL OBSERVATORY VECTOR ASTROMETRY SOFTWARE               *
*                                                                      *
*                            G. H. KAPLAN                              *
*                        U.S. NAVAL OBSERVATORY                        *
*                                                                      *
************************************************************************
 
 
 
      SUBROUTINE PLACE ( TJD, OBJECT, LOCATN, ICOORD, STAR, OBSERV,
     .                   SKYPOS )
*
*     THIS SUBROUTINE COMPUTES THE APPARENT DIRECTION OF A STAR OR SOLAR
*     SYSTEM BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE
*     SYSTEM.  BASED ON KAPLAN, ET AL. (1989), ASTRONOMICAL JOURNAL 97,
*     1197-1210, WITH SOME ENHANCEMENTS FROM KLIONER (2003),
*     ASTRONOMICAL JOURNAL 125, 1580-1597.
*
*          TJD    = TT JULIAN DATE FOR PLACE (IN)
*          OBJECT = CHARACTER STRING IDENTIFYING OBJECT OF INTEREST (IN)
*                   FOR SOLAR SYSTEM
*                   BODY,             SPECIFY THE NAME USING ALL UPPER-
*                                     CASE LETTERS ('SUN', 'MOON',
*                                     'JUPITER', ETC.),
*                                     - OR -
*                                     SPECIFY THE BODY ID NUMBER
*                                     IN A 4-CHARACTER STRING OF THE
*                                     FORM '=NNN', WHERE NNN IS THE
*                                     BODY ID NUMBER
*                   FOR STAR,         PROVIDE A BLANK STRING, THE WORD
*                                     'STAR', OR ANY STRING BEGINNING
*                                     WITH '*'
*          LOCATN = INTEGER CODE SPECIFYING LOCATION OF OBSERVER (IN)
*                   SET LOCATN=0 FOR OBSERVER AT GEOCENTER
*                   SET LOCATN=1 FOR OBSERVER ON SURFACE OF EARTH
*                   SET LOCATN=2 FOR OBSERVER ON NEAR-EARTH SPACECRAFT
*          ICOORD = INTEGER CODE SPECIFYING COORDINATE SYSTEM OF OUTPUT
*                   POSITION (IN)
*                   SET ICOORD=0 FOR GCRS (OR 'LOCAL GCRS')
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=2 FOR TRUE EQUATOR AND CIO OF DATE
*                   SET ICOORD=3 FOR ASTROMETRIC COORDINATES, I.E.,
*                                WITHOUT LIGHT DEFLECTION OR ABERRATION
*          STAR   = ARRAY OF CATALOG DATA FOR STAR (IN)
*                   (NOT USED IF SOLAR SYSTEM BODY REQUESTED)
*                   STAR(1) = ICRS RIGHT ASCENSION IN HOURS
*                   STAR(2) = ICRS DECLINATION IN DEGREES
*                   STAR(3) = ICRS PROPER MOTION IN RA IN
*                             MILLIARCSECONDS/YEAR
*                   STAR(4) = ICRS PROPER MOTION IN DEC IN
*                             MILLIARCSECONDS/YEAR
*                   STAR(5) = PARALLAX IN MILLIARCSECONDS
*                   STAR(6) = RADIAL VELOCITY IN KILOMETERS/SECOND
*                   FURTHER STAR ARRAY ELEMENTS ARE NOT USED HERE
*                   BUT ARE RESERVED FOR FUTURE USE
*          OBSERV = ARRAY OF DATA SPECIFYING LOCATION OF OBSERVER (IN)
*                   (NOT USED IF LOCATN=0)
*                   FOR LOCATN=1,
*                   OBSERV(1) = GEODETIC LONGITUDE (WGS-84) OF OBSERVER
*                               (EAST +) IN DEGREES
*                   OBSERV(2) = GEODETIC LATITUDE (WGS-84) OF OBSERVER
*                               (NORTH +) IN DEGREES
*                   OBSERV(3) = HEIGHT OF OBSERVER ABOVE ELLIPSOID
*                               IN METERS
*                   OBSERV(4) = VALUE OF DELTA-T IN SECONDS
*                               (DELTA-T=TT-UT1)
*                   OBSERV(5) = (NOT USED, RESERVED FOR FUTURE USE)
*                   OBSERV(6) = (NOT USED, RESERVED FOR FUTURE USE)
*                   FOR LOCATN=2,
*                   OBSERV(1) = GEOCENTRIC X IN KILOMETERS
*                   OBSERV(2) = GEOCENTRIC Y IN KILOMETERS
*                   OBSERV(3) = GEOCENTRIC Z IN KILOMETERS
*                   OBSERV(4) = GEOCENTRIC X-DOT IN KILOMETERS/SECOND
*                   OBSERV(5) = GEOCENTRIC Y-DOT IN KILOMETERS/SECOND
*                   OBSERV(6) = GEOCENTRIC Z-DOT IN KILOMETERS/SECOND
*                   WITH RESPECT TO TRUE EQUATOR AND EQUINOX OF DATE
*          SKYPOS = ARRAY OF OUTPUT DATA SPECIFYING OBJECT'S PLACE
*                   ON THE SKY AT TIME TJD, WITH RESPECT TO THE
*                   SPECIFIED OUTPUT COORDINATE SYSTEM (OUT)
*                   SKYPOS(1) = X, DIMENSIONLESS      UNIT VECTOR
*                   SKYPOS(2) = Y, DIMENSIONLESS      TOWARD OBJECT
*                   SKYPOS(3) = Z, DIMENSIONLESS
*                   SKYPOS(4) = APPARENT, TOPOCENTRIC, OR ASTROMETRIC
*                               RIGHT ASCENSION IN HOURS
*                   SKYPOS(5) = APPARENT, TOPOCENTRIC, OR ASTROMETRIC
*                               DECLINATION IN DEGREES
*                   SKYPOS(6) = TRUE (GEOMETRIC, EUCLIDIAN) DISTANCE
*                               TO SOLAR SYSTEM BODY IN AU AT TIME TJD,
*                               OR 0.D0 FOR STAR
*                   SKYPOS(7) = RADIAL VELOCITY IN KILOMETERS/SECOND
*                   FURTHER SKYPOS ARRAY ELEMENTS ARE NOT USED HERE
*                   BUT ARE RESERVED FOR FUTURE USE
*
*     NOTE 1: VALUES OF LOCATN AND ICOORD FOR VARIOUS STANDARD KINDS
*     OF PLACE:
*     LOCATN=0 AND ICOORD=1 APPARENT PLACE
*     LOCATN=1 AND ICOORD=1 TOPOCENTRIC PLACE
*     LOCATN=0 AND ICOORD=0 VIRTUAL PLACE
*     LOCATN=1 AND ICOORD=0 LOCAL PLACE
*     LOCATN=0 AND ICOORD=3 ASTROMETRIC PLACE
*     LOCATN=1 AND ICOORD=3 TOPOCENTRIC ASTROMETRIC PLACE
*
*     NOTE 2: ARRAYS STAR AND SKYPOS MAY BE EXPANDED IN THE FUTURE, AND
*     THIS CAN BE ALLOWED FOR IN THE CALLING CODE BY DIMENSIONING
*     THESE ARRAYS WITH 20 AND 10 ELEMENTS, RESPECTIVELY, EVEN THOUGH
*     ELEMENTS BEYOND STAR(6) AND SKYPOS(7) ARE NOT NOW REFERRED TO IN
*     THIS SUBROUTINE.
*
*     NOTE 3: IF LOCATN=1 AND OBSERV(4)=0.D0, THE VALUE OF DELTA-T WILL
*     BE OBTAINED FROM GETDT, WHICH PROVIDES THE LAST VALUE OF DELTA-T
*     DEFINED BY THE USER VIA CALL TO SETDT.
*
*     NOTE 4: SKYPOS(7), THE RADIAL VELOCITY, IS THE PREDICTED
*     RADIAL VELOCITY MEASURE (Z) TIMES THE SPEED OF LIGHT, AN
*     INHERENTLY SPECTROSCOPIC MEASURE.  FOR A STAR, IT
*     INCLUDES ALL EFFECTS, SUCH AS GRAVITATIONAL RED SHIFT,
*     CONTAINED IN THE CATALOG BARYCENTRIC RADIAL VELOCITY MEASURE,
*     WHICH IS ASSUMED GIVEN IN STAR(6).  FOR A SOLAR SYSTEM
*     BODY, IT APPLIES TO A FICTITIOUS EMITTER AT THE CENTER OF THE
*     OBSERVED OBJECT, ASSUMED MASSLESS (NO GRAVITATIONAL RED SHIFT),
*     AND DOES NOT IN GENERAL APPLY TO REFLECTED LIGHT.
*
*
 
* --- INITIAL DECLARATIONS---------------------------------------------
 
      IMPLICIT NONE
      INTEGER LOCATN,ICOORD,NTIMES,IEARTH,ISUN,IDBODY,IERR,LOC,J,KCIO,
     .     IDSS
      DOUBLE PRECISION TJD,STAR,OBSERV,SKYPOS,
     .     T0,TLAST1,TLAST2,TTJD,TDBJD,C,X,SECDIF,TLIGHT,DIS,DT,
     .     FRLIMB,RCIO,PEB,VEB,PSB,VSB,POG,VOG,POB,VOB,
     .     POS1,VEL1,POS2,POS3,POS4,POS5,POS6,POS7,POS8,
     .     PX,PY,PZ,RA,DEC,RVS,RVD,RV,DABS,DSQRT
      CHARACTER*(*) OBJECT
 
      DIMENSION STAR(*), OBSERV(6), SKYPOS(*),
     .     PEB(3), VEB(3), PSB(3), VSB(3),
     .     POG(3), VOG(3), POB(3), VOB(3),
     .     POS1(3), VEL1(3), POS2(3), POS3(3), POS4(3), POS5(3),
     .     POS6(3), POS7(3), POS8(3),
     .     PX(3), PY(3), PZ(3), RVS(3), RVD(3)
 
      SAVE
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST1, TLAST2 / 0.D0, 0.D0 /,   NTIMES / 0 /
 
   3  FORMAT ( ' PLACE: CANNOT OBTAIN COORDINATES OF ', A, ' AT JD ',
     .     F10.1 )
   4  FORMAT ( ' PLACE: WILL NOT PROCESS EARTH AS OBSERVED OBJECT ',
     .     'EXCEPT WHEN LOCATN=2' )
 
* --- GET CONSTANTS, FIRST TIME ONLY ----------------------------------
 
      NTIMES = NTIMES + 1
 
      IF ( NTIMES .EQ. 1 ) THEN
          IEARTH = IDSS ( 'EARTH' )
          ISUN = IDSS ( 'SUN' )
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ( 'C(AU/DAY)', 1.D0,   C )
      END IF
 
* --- CHECK ON EARTH AS AN OBSERVED OBJECT ----------------------------
 
      IF ( OBJECT .EQ. 'EARTH' .AND. LOCATN .NE. 2 ) THEN
          WRITE ( *, 4 )
          GO TO 70
      END IF
 
* --- GET POSITION AND VELOCITY OF EARTH (GEOCENTER) AND SUN ----------
 
      IF ( DABS ( TJD - TLAST1 ) .GT. 1.D-8 ) THEN
 
*         COMPUTE TDBJD, THE TDB JULIAN DATE CORRESPONDING TO TTJD
          TTJD = TJD
          TDBJD = TJD
          CALL TIMES ( TDBJD, X,   SECDIF )
          TDBJD = TTJD + SECDIF / 86400.D0
 
*         GET POSITION AND VELOCITY OF THE EARTH WRT BARYCENTER OF
*         SOLAR SYSTEM, IN ICRS
          CALL SOLSYS ( TDBJD, IEARTH, 0,   PEB, VEB, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) 'EARTH', TJD
              GO TO 70
          END IF
 
*         GET POSITION AND VELOCITY OF THE SUN WRT BARYCENTER OF
*         SOLAR SYSTEM, IN ICRS
          CALL SOLSYS ( TDBJD, ISUN, 0,   PSB, VSB, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) 'SUN', TJD
              GO TO 70
          END IF
 
          TLAST1 = TJD
 
      END IF
 
* --- GET POSITION AND VELOCITY OF OBSERVER ---------------------------
 
      IF ( LOCATN .EQ. 1 .OR. LOCATN .EQ. 2 ) THEN
 
*         FOR TOPOCENTRIC PLACE, GET GEOCENTRIC POSITION AND VELOCITY
*         VECTORS OF OBSERVER
          CALL GEOPOS ( TTJD, LOCATN, OBSERV,   POG, VOG )
          LOC = 1
 
      ELSE
 
*         FOR GEOCENTRIC PLACE, THERE IS NOTHING TO DO
          DO 25 J = 1, 3
             POG(J) = 0.D0
             VOG(J) = 0.D0
  25      CONTINUE
          LOC = 0
 
      END IF
 
*     COMPUTE POSITION AND VELOCITY OF OBSERVER WRT BARYCENTER OF
*     SOLAR SYSTEM (GALILEAN TRANSFORMATION FROM GCRS to BCRS)
      DO 30 J=1,3
          POB(J) = PEB(J) + POG(J)
          VOB(J) = VEB(J) + VOG(J)
  30  CONTINUE
 
* --- FIND GEOMETRIC POSITION OF OBSERVED OBJECT ----------------------
 
      IF ( OBJECT .EQ. 'STAR' .OR. OBJECT .EQ. ' ' .OR.
     .     OBJECT(1:1) .EQ. '*' ) THEN
 
*         OBSERVED OBJECT IS STAR
          IDBODY = -9999
 
*         GET POSITION OF STAR UPDATED FOR ITS SPACE MOTION
          CALL VECTRS ( STAR(1), STAR(2), STAR(3), STAR(4),
     .                  STAR(5), STAR(6),   POS1, VEL1 )
          CALL DLIGHT ( POS1, POB,   DT )
          CALL PROPMO ( T0, POS1, VEL1, TDBJD + DT,   POS2 )
*         GET POSITION OF STAR WRT OBSERVER (CORRECTED FOR PARALLAX)
          CALL GEOCEN ( POS2, POB,   POS3, TLIGHT )
          DIS = 0.D0
 
      ELSE
 
*         OBSERVED OBJECT IS SOLAR SYSTEM BODY
 
*         GET ID NUMBER OF BODY
          IF ( OBJECT(1:1) .EQ. '=' ) THEN
              READ ( OBJECT, '(1X,I3)' ) IDBODY
          ELSE
              IDBODY = IDSS ( OBJECT )
              IF ( IDBODY .EQ. -9999 ) THEN
                  WRITE ( *, 3 ) OBJECT, TJD
                  GO TO 70
              END IF
          END IF
 
*         GET POSITION OF BODY WRT BARYCENTER OF SOLAR SYSTEM
          CALL SOLSYS ( TDBJD, IDBODY, 0,   POS1, VEL1, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) OBJECT, TJD
              GO TO 70
          END IF
 
*         GET POSITION OF BODY WRT OBSERVER, AND TRUE (EUCLIDIAN)
*         DISTANCE
          CALL GEOCEN ( POS1, POB,   POS2, TLIGHT )
          DIS = TLIGHT * C
 
*         GET POSITION OF BODY WRT OBSERVER, ANTEDATED FOR LIGHT-TIME
          CALL LITTIM ( TDBJD, IDBODY, POB, TLIGHT,   POS3, TLIGHT )
 
      END IF
 
* --- APPLY GRAVITATIONAL DEFLECTION OF LIGHT AND ABERRATION ----------
 
      IF ( ICOORD .EQ. 3 ) THEN
 
*         THESE CALCULATIONS ARE SKIPPED FOR ASTROMETRIC PLACE
          DO 40 J = 1, 3
              POS5(J) = POS3(J)
  40      CONTINUE
 
      ELSE
 
*         VARIABLE LOC DETERMINES WHETHER EARTH DEFLECTION IS INCLUDED
          IF ( LOC .EQ. 1 ) THEN
              CALL LIMANG ( POS3, POG,   X, FRLIMB )
              IF ( FRLIMB .LT. 0.8D0 ) LOC = 0
          END IF
 
*         COMPUTE GRAVITATIONAL DEFLECTION AND ABERRATION
          CALL GRVDEF ( TDBJD, LOC, POS3, POB,   POS4 )
          CALL ABERAT ( POS4, VOB, TLIGHT,   POS5 )
*         POSITION VECTOR IS NOW IN GCRS
 
      END IF
 
* --- TRANSFORM, IF NECESSARY, TO OUTPUT COORDINATE SYSTEM ------------
 
      IF ( ICOORD .EQ. 1 ) THEN
 
*         TRANSFORM TO EQUATOR AND EQUINOX OF DATE
          CALL FRAME  ( POS5, 1,   POS6 )
          CALL PRECES ( T0, POS6, TDBJD,   POS7 )
          CALL NUTATE ( TDBJD, POS7,   POS8 )
 
      ELSE IF ( ICOORD .EQ. 2 ) THEN
 
*         TRANSFORM TO EQUATOR AND CIO OF DATE
          IF ( DABS ( TDBJD - TLAST2 ) .GT. 1.D-8 ) THEN
*             OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*             INTERMEDIATE SYSTEM
              CALL CIOLOC ( TDBJD,   RCIO, KCIO )
              CALL CIOBAS ( TDBJD, RCIO, KCIO,   PX, PY, PZ )
              TLAST2 = TDBJD
          END IF
*         TRANSFORM POSITION VECTOR TO CELESTIAL INTERMEDIATE SYSTEM
          POS8(1) = PX(1) * POS5(1) + PX(2) * POS5(2) + PX(3) * POS5(3)
          POS8(2) = PY(1) * POS5(1) + PY(2) * POS5(2) + PY(3) * POS5(3)
          POS8(3) = PZ(1) * POS5(1) + PZ(2) * POS5(2) + PZ(3) * POS5(3)
 
      ELSE
 
*         NO TRANSFORMATION -- KEEP COORDINATES IN GCRS
*         (OR ICRS FOR ASTROMETRIC PLACE)
          DO 50 J = 1, 3
              POS8(J) = POS5(J)
  50      CONTINUE
 
      END IF
 
* --- GET RADIAL VELOCITY ---------------------------------------------
 
*     SET UP STAR DATA, IF APPLICABLE
      IF ( IDBODY .EQ. -9999 ) THEN
          RVS(1) = STAR(1)
          RVS(2) = STAR(2)
          RVS(3) = STAR(6)
          IF ( STAR(5) .LE. 0.D0 ) THEN
              VEL1(1) = 0.D0
              VEL1(2) = 0.D0
              VEL1(3) = 0.D0
          END IF
      ELSE
          RVS(1) = 0.D0
          RVS(2) = 0.D0
          RVS(3) = 0.D0
      END IF
 
*     COMPUTE DISTANCES: OBSERVER-GEOCENTER, OBSERVER-SUN, OBJECT-SUN
      RVD(1) = DSQRT ( ( POB(1)  - PEB(1) ) ** 2
     .               + ( POB(2)  - PEB(2) ) ** 2
     .               + ( POB(3)  - PEB(3) ) ** 2 )
      RVD(2) = DSQRT ( ( POB(1)  - PSB(1) ) ** 2
     .               + ( POB(2)  - PSB(2) ) ** 2
     .               + ( POB(3)  - PSB(3) ) ** 2 )
      RVD(3) = DSQRT ( ( POS1(1) - PSB(1) ) ** 2
     .               + ( POS1(2) - PSB(2) ) ** 2
     .               + ( POS1(3) - PSB(3) ) ** 2 )
 
      CALL RADVL ( POS3, VEL1, VOB, RVS, RVD,   RV )
 
* --- FINISH UP -------------------------------------------------------
 
      CALL ANGLES ( POS8,   RA, DEC )
 
      X = DSQRT ( POS8(1)**2 + POS8(2)**2 + POS8(3)**2 )
 
      DO 60 J = 1, 3
          SKYPOS(J) = POS8(J) / X
  60  CONTINUE
 
      SKYPOS(4) = RA
      SKYPOS(5) = DEC
      SKYPOS(6) = DIS
      SKYPOS(7) = RV
 
      CALL SETVEC ( POS8 )
 
  70  RETURN
 
      END
 
 
 
      SUBROUTINE PLACES
*
*     THE ENTRIES TO THIS SUBROUTINE PROVIDE 'FRONT ENDS' TO
*     SUBROUTINE PLACE, TAILORED TO SPECIFIC PLACE TYPES.  THEY
*     PROVIDE COMPATIBILITY WITH PREVIOUSLY SUPPORTED CALLING
*     SEQUENCES.
*
*
      IMPLICIT NONE
      INTEGER L,N,LOCATN,ICOORD
      DOUBLE PRECISION TJD,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL,
     .     UJD,GLON,GLAT,HT,RA,DEC,DIS,
     .     TTJD,GAST,DELTAT,STAR,OBSERV,SKYPOS,DMOD
      CHARACTER*4 OBJECT
 
      DIMENSION STAR(20), OBSERV(6), SKYPOS(10)
 
      SAVE
 
      DATA STAR, OBSERV, SKYPOS / 36 * 0.D0 /
 
 
* --- GEOCENTRIC PLACES OF STARS --------------------------------------
*
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          TJD    = TT JULIAN DATE FOR APPARENT PLACE (IN)
*          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
*                   (NO LONGER USED)
*          RAI    = ICRS RIGHT ASCENSION IN HOURS (IN)
*          DECI   = ICRS DECLINATION IN DEGREES (IN)
*          PMRA   = ICRS PROPER MOTION IN RA IN MILLIARCSECONDS/YEAR
*                   (IN)
*          PMDEC  = ICRS PROPER MOTION IN DEC IN MILLIARCSECONDS/YEAR
*                   (IN)
*          PARLAX = PARALLAX IN MILLIARCSECONDS (IN)
*          RADVEL = RADIAL VELOCITY IN KILOMETERS/SECOND (IN)
*          RA     = APPARENT RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = APPARENT DECLINATION IN DEGREES (OUT)
*
*     NOTE: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS GCRS FOR
*     VPSTAR, ICRS FOR ASSTAR, AND EQUATOR AND EQUINOX OF DATE FOR
*     APSTAR.
*
*
*     THIS ENTRY PROVIDES THE APPARENT PLACE OF A STAR.
      ENTRY APSTAR (TJD,N,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL, RA,DEC)
      LOCATN = 0
      ICOORD = 1
      GO TO 10
 
*     THIS ENTRY PROVIDES THE VIRTUAL PLACE OF A STAR.
      ENTRY VPSTAR (TJD,N,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL, RA,DEC)
      LOCATN = 0
      ICOORD = 0
      GO TO 10
 
*     THIS ENTRY PROVIDES THE ASTROMETRIC PLACE OF A STAR.
      ENTRY ASSTAR (TJD,N,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL, RA,DEC)
      LOCATN = 0
      ICOORD = 3
 
  10  TTJD = TJD
      OBJECT = 'STAR'
      STAR(1) = RAI
      STAR(2) = DECI
      STAR(3) = PMRA
      STAR(4) = PMDEC
      STAR(5) = PARLAX
      STAR(6) = RADVEL
 
      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
 
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
 
      RETURN
 
 
* --- TOPOCENTRIC PLACES OF STARS -------------------------------------
*
*     EACH OF THESE ENTRIES CAN BE CALLED ONLY AFTER A CALL TO THE
*     CORRESPONDING GEOCENTRIC ENTRY, WHICH SUPPLIES SOME REQUIRED DATA.
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          UJD    = UT1 JULIAN DATE FOR TOPOCENTRIC PLACE (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = TOPOCENTRIC DECLINATION IN DEGREES (OUT)
*
*     NOTE 1: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS 'LOCAL GCRS'
*     FOR LPSTAR, AND EQUATOR AND EQUINOX OF DATE FOR TPSTAR.
*
*     NOTE 2: UJD CAN ALSO BE GREENWICH APPARENT SIDEREAL TIME IN HOURS,
*     EQUIVALENT TO UT1 JULIAN DATE, BUT THIS OPTION WILL NOT BE
*     SUPPORTED INDEFINITELY.  ADVISE USING UJD = UT1 JULIAN DATE ONLY.
*
*
*     THIS ENTRY PROVIDES THE TOPOCENTRIC PLACE OF A STAR.
      ENTRY TPSTAR (UJD,GLON,GLAT,HT, RA,DEC)
      LOCATN = 1
      ICOORD = 1
      GO TO 20
 
*     THIS ENTRY PROVIDES THE LOCAL PLACE OF A STAR.
      ENTRY LPSTAR (UJD,GLON,GLAT,HT, RA,DEC)
      LOCATN = 1
      ICOORD = 0
 
  20  IF ( UJD .GT. 100.D0 ) THEN
          DELTAT = ( TTJD - UJD ) * 86400.D0
      ELSE
          GAST = DMOD ( UJD, 24.D0 )
          IF ( GAST .LT. 0.D0 ) GAST = GAST + 24.D0
          CALL PLACST ( GAST )
          DELTAT = 0.D0
      END IF
      OBSERV(1) = GLON
      OBSERV(2) = GLAT
      OBSERV(3) = HT
      OBSERV(4) = DELTAT
 
      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
 
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
 
      RETURN
 
 
* --- GEOCENTRIC PLACES OF PLANETS ------------------------------------
*
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          TJD    = TT JULIAN DATE FOR APPARENT PLACE (IN)
*          L      = BODY IDENTIFICATION NUMBER FOR DESIRED PLANET (IN)
*          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
*                   (NO LONGER USED)
*          RA     = APPARENT RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = APPARENT DECLINATION IN DEGREES (OUT)
*          DIS    = TRUE DISTANCE FROM EARTH TO PLANET IN AU (OUT)
*
*     NOTE: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS GCRS FOR
*     VPPLAN, ICRS FOR ASPLAN, AND EQUATOR AND EQUINOX OF DATE FOR
*     APPLAN.
*
*     NOTE: 'PLANET' IS USED GENERICALLY FOR ANY SOLAR SYSTEM BODY.
*
*     THIS ENTRY PROVIDES THE APPARENT PLACE OF A PLANET.
      ENTRY APPLAN (TJD,L,N, RA,DEC,DIS)
      LOCATN = 0
      ICOORD = 1
      GO TO 30
 
*     THIS ENTRY PROVIDES THE VIRTUAL PLACE OF A PLANET.
      ENTRY VPPLAN (TJD,L,N, RA,DEC,DIS)
      LOCATN = 0
      ICOORD = 0
      GO TO 30
 
*     THIS ENTRY PROVIDES THE ASTROMETRIC PLACE OF A PLANET.
      ENTRY ASPLAN (TJD,L,N, RA,DEC,DIS)
      LOCATN = 0
      ICOORD = 3
 
  30  TTJD = TJD
      WRITE ( OBJECT, '(''='',I3)' ) L
 
      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
 
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
      DIS = SKYPOS(6)
 
      RETURN
 
 
* --- TOPOCENTRIC PLACES OF PLANETS -----------------------------------
*
*     EACH OF THESE ENTRIES CAN BE CALLED ONLY AFTER A CALL TO THE
*     CORRESPONDING GEOCENTRIC ENTRY, WHICH SUPPLIES SOME REQUIRED DATA.
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          UJD    = UT1 JULIAN DATE FOR TOPOCENTRIC PLACE (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = TOPOCENTRIC DECLINATION IN DEGREES (OUT)
*          DIS    = TRUE DISTANCE FROM OBSERVER TO PLANET IN AU (OUT)
*
*     NOTE 1: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS 'LOCAL GCRS'
*     FOR LPPLAN, AND EQUATOR AND EQUINOX OF DATE FOR TPPLAN.
*
*     NOTE 2: UJD CAN ALSO BE GREENWICH APPARENT SIDEREAL TIME IN HOURS,
*     EQUIVALENT TO UT1 JULIAN DATE, BUT THIS OPTION WILL NOT BE
*     SUPPORTED INDEFINITELY.  ADVISE USING UJD = UT1 JULIAN DATE ONLY.
*
*
*     THIS ENTRY PROVIDES THE TOPOCENTRIC PLACE OF A PLANET.
      ENTRY TPPLAN (UJD,GLON,GLAT,HT, RA,DEC,DIS)
      LOCATN = 1
      ICOORD = 1
      GO TO 40
 
*     THIS ENTRY PROVIDES THE LOCAL PLACE OF A PLANET.
      ENTRY LPPLAN (UJD,GLON,GLAT,HT, RA,DEC,DIS)
      LOCATN = 1
      ICOORD = 0
 
  40  IF ( UJD .GT. 100.D0 ) THEN
          DELTAT = ( TTJD - UJD ) * 86400.D0
      ELSE
          GAST = DMOD ( UJD, 24.D0 )
          IF ( GAST .LT. 0.D0 ) GAST = GAST + 24.D0
          CALL PLACST ( GAST )
          DELTAT = 0.D0
      END IF
      OBSERV(1) = GLON
      OBSERV(2) = GLAT
      OBSERV(3) = HT
      OBSERV(4) = DELTAT
 
      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
 
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
      DIS = SKYPOS(6)
 
      RETURN
 
      END
 
 
 
      SUBROUTINE MPSTAR (TJD,N,RA,DEC, RAI,DECI)
*
*     THIS SUBROUTINE COMPUTES THE ICRS POSITION OF A STAR,
*     GIVEN ITS APPARENT PLACE AT DATE TJD.  PROPER MOTION, PARALLAX,
*     AND RADIAL VELOCITY ARE ASSUMED TO BE ZERO.
*
*          TJD    = TT JULIAN DATE OF APPARENT PLACE (IN)
*          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
*                   (NO LONGER USED)
*          RA     = APPARENT RIGHT ASCENSION IN HOURS, REFERRED TO
*                   TRUE EQUATOR AND EQUINOX OF DATE (IN)
*          DEC    = APPARENT DECLINATION IN DEGREES, REFERRED TO
*                   TRUE EQUATOR OF DATE (IN)
*          RAI    = ICRS RIGHT ASCENSION IN HOURS (OUT)
*          DECI   = ICRS DECLINATION IN DEGREES (OUT)
*
*
      DOUBLE PRECISION TJD,RA,DEC,RAI,DECI,T0,T1,RAINEW,DCINEW,
     .     RAIOLD,DCIOLD,STAR,OBSERV,SKYPOS,R,D,P,V,DELRA,DELDEC,
     .     DABS
      DIMENSION STAR(20), OBSERV(6), SKYPOS(10), P(3), V(3)
      SAVE
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA STAR, OBSERV, SKYPOS / 36 * 0.D0 /
 
   3  FORMAT ( ' MPSTAR: NO CONVERGENCE AT COORDINATES ',
     .     F9.5, 1X, SP,F9.4, ' AT JD ', SS, F10.1 )
 
      T1 = TJD
 
*     GET INITIAL APPROXIMATION
      ITER = 0
      CALL VECTRS (RA,DEC,0.D0,0.D0,0.D0,0.D0,P,V)
      CALL PRECES (T1,P,T0,V)
      CALL ANGLES (V,RAINEW,DCINEW)
 
*     ITERATIVELY FIND ICRS COORDINATES THAT PRODUCE INPUT
*     APPARENT PLACE OF STAR AT DATE TJD
   20 ITER = ITER + 1
      RAIOLD = RAINEW
      DCIOLD = DCINEW
      STAR(1) = RAIOLD
      STAR(2) = DCIOLD
      STAR(3) = 0.D0
      STAR(4) = 0.D0
      STAR(5) = 0.D0
      STAR(6) = 0.D0
      CALL PLACE (T1,'STAR',0,1,STAR,OBSERV, SKYPOS)
      R = SKYPOS(4)
      D = SKYPOS(5)
      DELRA = R - RA
      DELDEC = D - DEC
      IF (DELRA.LT.-12.D0) DELRA = DELRA + 24.D0
      IF (DELRA.GT.+12.D0) DELRA = DELRA - 24.D0
      RAINEW = RAIOLD - DELRA
      DCINEW = DCIOLD - DELDEC
      IF (ITER.GT.30) THEN
          WRITE ( *, 3 ) RA, DEC, TJD
          GO TO 40
      END IF
      IF (DABS(RAINEW-RAIOLD).GT.1.D-12) GO TO 20
      IF (DABS(DCINEW-DCIOLD).GT.1.D-11) GO TO 20
 
   40 RAI = RAINEW
      DECI = DCINEW
      IF (RAI.LT. 0.D0) RAI = RAI + 24.D0
      IF (RAI.GE.24.D0) RAI = RAI - 24.D0
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL VECTRS (RAI,DECI,0.D0,0.D0,0.D0,0.D0,P,V)
      CALL SETVEC (P)
 
   50 RETURN
 
      END
 
 
 
      SUBROUTINE SIDTIM ( TJDH, TJDL, K,   GST )
*
*     THIS SUBROUTINE COMPUTES THE GREENWICH SIDEREAL TIME
*     (EITHER MEAN OR APPARENT) AT JULIAN DATE TJDH + TJDL.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          K      = TIME SELECTION CODE (IN)
*                   SET K=0 FOR GREENWICH MEAN SIDEREAL TIME
*                   SET K=1 FOR GREENWICH APPARENT SIDEREAL TIME
*          GST    = GREENWICH (MEAN OR APPARENT) SIDEREAL TIME
*                   IN HOURS (OUT)
*
*     NOTE:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*
      DOUBLE PRECISION TJDH,TJDL,GST,PI,DEGCON,DELTAT,
     .     T0,UTJD,TTJD,TDBJD,SECDIF,A,THETA,RCIO,
     .     UNITX,W1,W2,X,Y,Z,EQ,HAEQ,EE,DMOD,DATAN2
      DIMENSION UNITX(3), W1(3), W2(3), X(3), Y(3), Z(3), EQ(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGCON = 180.D0 / PI           )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA UNITX / 1.D0, 0.D0, 0.D0 /
 
   3  FORMAT ( ' SIDTIM ERROR: CANNOT RETURN SIDEREAL TIME FOR ',
     .     'JD ', F10.1 )
 
      CALL GETDT ( DELTAT )
 
*     TIME ARGUMENT FOR PRECESSION AND NUTATION COMPONENTS OF SIDEREAL
*     TIME IS TDB
      UTJD = TJDH + TJDL
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, A,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0
 
*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )
 
      IF ( MODE .GE. 2 ) THEN
*          EQUINOX-BASED MODE
*          SEE USNO CIRCULAR 179, SECTION 2.6.2
 
*         GET -1 TIMES THE MEAN OR TRUE RIGHT ASCENSION OF THE CIO
          CALL EQXRA ( TDBJD, K,   RCIO )
*         GET EARTH ROTATION ANGLE
          CALL EROT ( TJDH, TJDL,   THETA )
*         COMBINE TO OBTAIN SIDEREAL TIME
          GST = DMOD ( THETA / 15.D0 - RCIO, 24.D0 )
          IF ( GST .LT. 0.D0 ) GST = GST + 24.D0
 
      ELSE
*         CIO-BASED MODE
*         SEE USNO CIRCULAR 179, SECTION 6.5.4
 
*         GET EARTH ROTATION ANGLE
          CALL EROT ( TJDH, TJDL,   THETA )
*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          IF ( RCIO .EQ. 99.D0 ) THEN
              WRITE ( *, 3 ) TDBJD
              GST = 99.D0
              GO TO 50
          END IF
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )
*         COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
          CALL NUTATE ( -TDBJD, UNITX,   W1 )
          CALL PRECES ( TDBJD, W1, T0,   W2 )
          CALL FRAME ( W2, -1,    EQ )
*         COMPUTE THE HOUR ANGLE OF THE EQUINOX WRT THE TIO MERIDIAN
*         (NEAR GREENWICH, BUT PASSES THROUGH THE CIP AND TIO)
          HAEQ = THETA - DATAN2 ( EQ(1)*Y(1) + EQ(2)*Y(2) + EQ(3)*Y(3),
     .                            EQ(1)*X(1) + EQ(2)*X(2) + EQ(3)*X(3) )
     .                   * DEGCON
 
*         FOR MEAN SIDEREAL TIME, OBTAIN THE EQUATION OF THE EQUINOXES
*         AND SUBTRACT IT
          IF ( K .EQ. 0 ) THEN
              CALL ETILT ( TDBJD,   A, A, EE, A, A )
              HAEQ = HAEQ - EE / 240.D0
          END IF
 
          HAEQ = DMOD ( HAEQ, 360.D0 ) / 15.D0
          IF ( HAEQ .LT. 0.D0 ) HAEQ = HAEQ + 24.D0
          GST = HAEQ
 
      END IF
 
  50  RETURN
 
      END
 
 
 
      SUBROUTINE CIORA ( TJD,   RACIO )
*
*     THIS SUBROUTINE COMPUTES THE TRUE RIGHT ASCENSION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO) AT A GIVEN TT JULIAN DATE.  THIS IS
*     -(EQUATION OF THE ORIGINS).
*
*          TJD    = TT JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, WITH RESPECT TO THE
*                   TRUE EQUINOX OF DATE, IN HOURS (+ OR -) (OUT)
*
*
      DOUBLE PRECISION TJD,RACIO,PI,DEGCON,T0,T,SECDIF,TDBJD,RCIO,
     .     UNITX,W1,W2,X,Y,Z,EQ,AZ,DATAN2
      DIMENSION UNITX(3), W1(3), W2(3), X(3), Y(3), Z(3), EQ(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGCON = 180.D0 / PI           )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA UNITX / 1.D0, 0.D0, 0.D0 /
 
   3  FORMAT ( ' CIORA ERROR: CANNOT RETURN CIO RA VALUE FOR JD ',
     .     F10.1 )
 
*     TDBJD IS THE TDB JULIAN DATE
      TDBJD = TJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TJD + SECDIF / 86400.D0
 
*     OBTAIN THE BASIS VECTORS, IN THE GCRS, FOR THE CELESTIAL
*     INTERMEDIATE SYSTEM DEFINED BY THE CIP (IN THE Z DIRECTION) AND
*     THE CIO (IN THE X DIRECTION)
      CALL CIOLOC ( TDBJD,   RCIO, KCIO )
      IF ( RCIO .EQ. 99.D0 ) THEN
          WRITE ( *, 3 ) TDBJD
          RACIO = 99.D0
          GO TO 50
      END IF
      CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )
 
*     COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
      CALL NUTATE ( -TDBJD, UNITX,   W1 )
      CALL PRECES ( TDBJD, W1, T0,   W2 )
      CALL FRAME ( W2, -1,    EQ )
 
*     COMPUTE THE INTERMEDIATE RA OF THE TRUE EQUINOX (EQUATION OF
*     THE ORIGINS)
      AZ = DATAN2 ( EQ(1)*Y(1) + EQ(2)*Y(2) + EQ(3)*Y(3),
     .              EQ(1)*X(1) + EQ(2)*X(2) + EQ(3)*X(3) ) * DEGCON
 
*     THE RA OF THE CIO IS MINUS THIS COORDINATE
      RACIO = -AZ / 15.D0
 
  50  RETURN
 
      END
 
 
 
      SUBROUTINE TERCEL ( TJDH, TJDL, XP, YP, VEC1,   VEC2 )
*
*     THIS SUBROUTINE ROTATES A VECTOR FROM THE TERRESTRIAL TO THE
*     CELESTIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
*     ITRS (A ROTATING EARTH-FIXED SYSTEM) TO THE GCRS (A LOCAL
*     SPACE-FIXED SYSTEM) BY APPLYING ROTATIONS FOR POLAR MOTION,
*     EARTH ROTATION, NUTATION, PRECESSION, AND THE DYNAMICAL-TO-GCRS
*     FRAME TIE.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
*                   SYSTEM) (IN)
*          VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
*                   SYSTEM) (OUT)
*
*     NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
*
*     NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*     NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
*     (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
*     SPECIAL OPTION WHERE THE OUTPUT VECTOR IS PRODUCED WITH RESPECT
*     TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH THE
*     TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
*     WORKS ONLY IN 'EQUINOX' MODE.
*
*     NOTE 4: INPUT PARAMETERS XP, YP WERE XPOLE, YPOLE IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION TJDH,TJDL,XP,YP,VEC1,VEC2,
     .     T0,DELTAT,UTJDH,UTJDL,UTJD,TTJD,TDBJD,T,SECDIF,
     .     GAST,RCIO,THETA,V1,V2,V3,V4,X,Y,Z
      DIMENSION VEC1(3), VEC2(3), V1(3), V2(3), V3(3), V4(3),
     .     X(3), Y(3), Z(3)
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
 
      CALL GETDT ( DELTAT )
 
      IF ( TJDH .GE. 0.D0 ) THEN
          UTJDH = TJDH
          UTJDL = TJDL
      ELSE
          UTJDH = TJDL
          UTJDL = 0.D0
      END IF
      UTJD = UTJDH + UTJDL
 
*     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0
 
*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )
 
      IF ( MODE .GE. 2 ) THEN
*         'EQUINOX' MODE
 
*         APPLY POLAR MOTION
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              V1(1) = VEC1(1)
              V1(2) = VEC1(2)
              V1(3) = VEC1(3)
          ELSE
              CALL WOBBLE ( TDBJD, XP, YP, VEC1,   V1 )
          END IF
 
*         APPLY EARTH ROTATION
          CALL SIDTIM ( UTJDH, UTJDL, 1,   GAST )
          CALL SPIN ( -GAST * 15.D0, V1,   V2 )
 
*         SPECIAL OPTION SKIPS REMAINING TRANSFORMATIONS
          IF ( TJDH .LT. 0.D0 ) THEN
              VEC2(1) = V2(1)
              VEC2(2) = V2(2)
              VEC2(3) = V2(3)
          ELSE
 
*         APPLY NUTATION AND PRECESSION
          CALL NUTATE ( -TDBJD, V2,   V3 )
          CALL PRECES ( TDBJD, V3, T0,   V4 )
 
*         APPLY FRAME-TIE MATRIX
          CALL FRAME ( V4, -1, VEC2 )
 
          END IF
 
      ELSE
*         'CIO-TIO-THETA' MODE
*         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
*         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
*         EQ. (3) AND (4).
 
*         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE TERRESTRIAL
*         INTERMEDIATE SYSTEM
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              V1(1) = VEC1(1)
              V1(2) = VEC1(2)
              V1(3) = VEC1(3)
          ELSE
              CALL WOBBLE ( TDBJD, XP, YP, VEC1,   V1 )
          END IF
 
*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )
 
*         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
*         THE VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
          CALL EROT ( UTJDH, UTJDL,   THETA )
          CALL SPIN ( -THETA, V1,   V2 )
 
*         TRANSFORM THE VECTOR FROM THE CELESTIAL INTERMEDIATE SYSTEM
*         TO THE GCRS
          VEC2(1) = X(1) * V2(1) + Y(1) * V2(2) + Z(1) * V2(3)
          VEC2(2) = X(2) * V2(1) + Y(2) * V2(2) + Z(2) * V2(3)
          VEC2(3) = X(3) * V2(1) + Y(3) * V2(2) + Z(3) * V2(3)
 
      END IF
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
  50  CALL SETVEC ( VEC2 )
 
      RETURN
 
      END
 
 
      SUBROUTINE CELTER ( TJDH, TJDL, XP, YP, VEC1,   VEC2 )
*
*     THIS SUBROUTINE ROTATES A VECTOR FROM THE CELESTIAL TO THE
*     TERRESTRIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
*     GCRS (A LOCAL SPACE-FIXED SYSTEM) TO THE ITRS (A ROTATING
*     EARTH-FIXED SYSTEM) BY APPLYING ROTATIONS FOR THE GCRS-TO-
*     DYNAMICAL FRAME TIE, PRECESSION, NUTATION, EARTH ROTATION,
*     AND POLAR MOTION.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
*                   SYSTEM) (IN)
*          VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
*                   SYSTEM) (OUT)
*
*     NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
*
*     NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*     NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
*     (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
*     SPECIAL OPTION WHERE THE INPUT VECTOR IS ASSUMED TO BE WITH
*     RESPECT TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH
*     THE TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
*     WORKS ONLY IN 'EQUINOX' MODE.
*
*
      DOUBLE PRECISION TJDH,TJDL,XP,YP,VEC1,VEC2,
     .     T0,DELTAT,UTJDH,UTJDL,UTJD,TTJD,TDBJD,T,SECDIF,
     .     GAST,RCIO,THETA,V1,V2,V3,V4,X,Y,Z
      DIMENSION VEC1(3), VEC2(3), V1(3), V2(3), V3(3), V4(3),
     .     X(3), Y(3), Z(3)
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
 
      CALL GETDT ( DELTAT )
 
      IF ( TJDH .GE. 0.D0 ) THEN
          UTJDH = TJDH
          UTJDL = TJDL
      ELSE
          UTJDH = TJDL
          UTJDL = 0.D0
      END IF
      UTJD = UTJDH + UTJDL
 
*     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0
 
*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )
 
      IF ( MODE .GE. 2 ) THEN
*         'EQUINOX' MODE
 
*         SPECIAL OPTION SKIPS INITIAL TRANSFORMATIONS
          IF ( TJDH .LT. 0.D0 ) THEN
              V3(1) = VEC1(1)
              V3(2) = VEC1(2)
              V3(3) = VEC1(3)
          ELSE
 
*         APPLY FRAME-TIE MATRIX
          CALL FRAME ( VEC1, 1, V1 )
 
*         APPLY PRECESSION AND NUTATION
          CALL PRECES ( T0, V1, TDBJD,   V2 )
          CALL NUTATE ( TDBJD, V2,   V3 )
 
          END IF
 
*         APPLY EARTH ROTATION
          CALL SIDTIM ( UTJDH, UTJDL, 1,   GAST )
          CALL SPIN ( GAST * 15.D0, V3,   V4 )
 
*         APPLY POLAR MOTION
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              VEC2(1) = V4(1)
              VEC2(2) = V4(2)
              VEC2(3) = V4(3)
          ELSE
              CALL WOBBLE ( -TDBJD, XP, YP, V4,   VEC2 )
          END IF
 
      ELSE
*         'CIO-TIO-THETA' MODE
*         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
*         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
*         EQ. (3) AND (4).
 
*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )
 
*         TRANSFORM THE VECTOR FROM THE GCRS TO THE
*         CELESTIAL INTERMEDIATE SYSTEM
          V1(1) = X(1) * VEC1(1) + X(2) * VEC1(2) + X(3) * VEC1(3)
          V1(2) = Y(1) * VEC1(1) + Y(2) * VEC1(2) + Y(3) * VEC1(3)
          V1(3) = Z(1) * VEC1(1) + Z(2) * VEC1(2) + Z(3) * VEC1(3)
 
*         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
*         THE VECTOR TO THE TERRESTRIAL INTERMEDIATE SYSTEM
          CALL EROT ( UTJDH, UTJDL,   THETA )
          CALL SPIN ( THETA, V1,   V2 )
 
*         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE ITRS
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              VEC2(1) = V2(1)
              VEC2(2) = V2(2)
              VEC2(3) = V2(3)
          ELSE
              CALL WOBBLE ( -TDBJD, XP, YP, V2,   VEC2 )
          END IF
 
      END IF
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
  50  CALL SETVEC ( VEC2 )
 
      RETURN
 
      END
 
 
 
 
      SUBROUTINE GETHIP ( RAH, DECH, PMRAH, PMDECH, PARXH, RVH,
     .                    RA2, DEC2, PMRA2, PMDEC2, PARX2, RV2 )
*
*     THIS SUBROUTINE CONVERTS HIPPARCOS DATA AT EPOCH J1991.25
*     TO EPOCH J2000.0.  TO BE USED ONLY FOR HIPPARCOS OR TYCHO STARS
*     WITH LINEAR SPACE MOTION.  BOTH INPUT AND OUTPUT DATA IS IN THE
*     ICRS.
*
*          RAH    = HIPPARCOS RIGHT ASCENSION IN DEGREES (IN)
*          DECH   = HIPPARCOS DECLINATION IN DEGREES (IN)
*          PMRAH  = HIPPARCOS PROPER MOTION IN RA
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PMDECH = HIPPARCOS PROPER MOTION IN DEC
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PARXH  = HIPPARCOS PARALLAX IN MILLIARCSECONDS (IN)
*          RVH    = RADIAL VELOCITY AT HIPPARCOS EPOCH
*                   IN KILOMETERS/SECOND (IN)
*          RA2    = RIGHT ASCENSION AT J2000.0 IN HOURS (OUT)
*          DEC2   = DECLINATION AT J2000.0 IN DEGREES (OUT)
*          PMRA2  = PROPER MOTION IN RA AT J2000.0
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PMDEC2 = PROPER MOTION IN DEC AT J2000.0
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PARX2  = PARALLAX AT J2000.0 IN MILLIARCSECONDS (OUT)
*          RV2    = RADIAL VELOCITY AT J2000.0 IN KILOMETERS/SECOND
*                   (OUT)
*
*     NOTE:  INPUT RA IS IN DEGREES, AS PER HIPPARCOS, BUT OUTPUT RA
*     IS IN HOURS.
*
*
      DOUBLE PRECISION RAH,DECH,PMRAH,PMDECH,PARXH,RVH,
     .     RA2,DEC2,PMRA2,PMDEC2,PARX2,RV2,EPOCH1,EPOCH2
 
      DATA EPOCH1, EPOCH2 / 2448349.0625D0, 2451545.0000D0 /
 
      CALL CATRAN ( 1,EPOCH1,RAH/15.D0,DECH,PMRAH,PMDECH,PARXH,RVH,
     .                EPOCH2,RA2,      DEC2,PMRA2,PMDEC2,PARX2,RV2 )
 
      RETURN
 
      END
 
 
 
      SUBROUTINE CATRAN ( IT,
     .                    DATE1, RA1, DEC1, PMRA1, PMDEC1, PARX1, RV1,
     .                    DATE2, RA2, DEC2, PMRA2, PMDEC2, PARX2, RV2 )
*
*     THIS SUBROUTINE TRANSFORMS A STAR'S CATALOG QUANTITIES FOR
*     A CHANGE OF EPOCH AND/OR EQUATOR AND EQUINOX.  IT CAN ALSO BE
*     USED TO ROTATE CATALOG QUANTITIES ON THE DYNAMICAL EQUATOR AND
*     EQUINOX OF J2000.0 TO THE ICRS OR VICE VERSA.
*
*          IT     = TRANSFORMATION OPTION (IN)
*                   SET IT=1 TO CHANGE EPOCH (SAME EQUATOR AND EQUINOX)
*                   SET IT=2 TO CHANGE EQUATOR AND EQUINOX (SAME EPOCH)
*                   SET IT=3 TO CHANGE EQUATOR AND EQUINOX AND EPOCH
*                   SET IT=4 TO CHANGE EQUATOR AND EQUINOX OF J2000.0
*                            TO ICRS
*                   SET IT=5 TO CHANGE ICRS TO EQUATOR AND EQUINOX OF
*                            J2000.0
*          DATE1  = TT JULIAN DATE, OR YEAR, OF ORIGINAL CATALOG
*                   DATA (THE FOLLOWING SIX ARGUMENTS) (IN)
*          RA1    = ORIGINAL MEAN RIGHT ASCENSION IN HOURS (IN)
*          DEC1   = ORIGINAL MEAN DECLINATION IN DEGREES (IN)
*          PMRA1  = ORIGINAL PROPER MOTION IN RA
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PMDEC1 = ORIGINAL PROPER MOTION IN DEC
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PARX1  = ORIGINAL PARALLAX IN MILLIARCSECONDS (IN)
*          RV1    = ORIGINAL RADIAL VELOCITY IN KILOMETERS/SECOND
*                   (IN)
*          DATE2  = TT JULIAN DATE, OR YEAR, FOR TRANSFORMED
*                   OUTPUT DATA (THE FOLLOWING SIX ARGUMENTS) (IN)
*          RA2    = TRANSFORMED MEAN RIGHT ASCENSION IN HOURS (OUT)
*          DEC2   = TRANSFORMED MEAN DECLINATION IN DEGREES (OUT)
*          PMRA2  = TRANSFORMED PROPER MOTION IN RA
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PMDEC2 = TRANSFORMED PROPER MOTION IN DEC
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PARX2  = TRANSFORMED PARALLAX IN MILLIARCSECONDS (OUT)
*          RV2    = TRANSFORMED RADIAL VELOCITY IN KILOMETERS/SECOND
*                   (OUT)
*
*     NOTE 1:  DATE1 AND DATE2 MAY BE SPECIFIED EITHER AS A JULIAN
*     DATE (E.G., 2433282.5D0) OR A JULIAN YEAR AND FRACTION
*     (E.G., 1950.0D0).  VALUES LESS THAN 10000 ARE ASSUMED TO
*     BE YEARS.  FOR IT=2 OR IT=3, EITHER DATE1 OR DATE2 MUST BE
*     2451545.0 OR 2000.0 (J2000.0).  FOR IT=4 AND IT=5, DATE1 AND
*     DATE2 ARE IGNORED.
*
*     NOTE 2:  IT=1 UPDATES THE STAR'S DATA TO ACCOUNT FOR
*     THE STAR'S SPACE MOTION BETWEEN THE FIRST AND SECOND DATES,
*     WITHIN A FIXED REFERENCE SYSTEM.  IT=2 APPLIES A ROTATION
*     OF THE REFERENCE SYSTEM CORRESPONDING TO PRECESSION BETWEEN
*     THE FIRST AND SECOND DATES, BUT LEAVES THE STAR FIXED IN SPACE.
*     IT=3 PROVIDES BOTH TRANSFORMATIONS.  IT=4 AND IT=5 PROVIDE A
*     A FIXED ROTATION ABOUT VERY SMALL ANGLES (<0.1 ARCSECOND) TO
*     TAKE DATA FROM THE DYNAMICAL SYSTEM OF J2000.0 TO THE ICRS (IT=4)
*     OR VICE VERSA (IT=5).
*
*     NOTE 3:  FOR IT=1, INPUT DATA CAN BE IN ANY FIXED REFERENCE
*     SYSTEM. FOR IT=2 OR IT=3, THIS SUBROUTINE ASSUMES THE INPUT DATA
*     IS IN THE DYNAMICAL SYSTEM AND PRODUCES OUTPUT IN THE DYNAMICAL
*     SYSTEM.  FOR IT=4, THE INPUT DATA MUST BE ON THE DYNAMICAL EQUATOR
*     AND EQUINOX OF J2000.0.  FOR IT=5, THE INPUT DATA MUST BE IN THE
*     ICRS.
*
*     NOTE 4:  THIS SUBROUTINE CANNOT BE PROPERLY USED TO BRING DATA
*     FROM OLD STAR CATALOGS INTO THE MODERN SYSTEM, BECAUSE
*     OLD CATALOGS WERE COMPILED USING A SET OF CONSTANTS THAT ARE
*     INCOMPATIBLE WITH MODERN VALUES.  IN PARTICULAR, IT SHOULD NOT
*     BE USED FOR CATALOGS WHOSE POSITIONS AND PROPER MOTIONS WERE
*     DERIVED BY ASSUMING A PRECESSION CONSTANT SIGNIFICANTLY DIFFERENT
*     FROM THE VALUE IMPLICIT IN SUBROUTINE PRECES.
*
*
      DOUBLE PRECISION DATE1,RA1,DEC1,PMRA1,PMDEC1,PARX1,RV1,
     .     DATE2,RA2,DEC2,PMRA2,PMDEC2,PARX2,RV2,
     .     PI,SECCON,AUKM,C,TJD1,POS1,VEL1,TJD2,POS2,VEL2,
     .     PARALX,DIST,R,D,CRA,SRA,CDC,SDC,K,PMR,PMD,RVL,
     .     XYPROJ,DCOS,DSIN,DATAN2
      INTEGER IT,J
 
      DIMENSION POS1(3), VEL1(3), POS2(3), VEL2(3)
 
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET LENGTH OF AU IN KILOMETERS
          CALL ASTCON ( 'AU', 1.D-3, AUKM )
*         GET C, THE SPEED OF LIGHT IN KILOMETERS/SECOND
          CALL ASTCON ( 'C', 1.D-3,   C )
      END IF
 
* --- IF NECESSARY, COMPUTE JULIAN DATES ------------------------------
 
*     SUBROUTINE USES TDB JULIAN DATES INTERNALLY, BUT NO
*     DISTINCTION BETWEEN TDB AND TT IS NECESSARY
 
      IF ( DATE1 .LT. 10000.D0 ) THEN
           TJD1 = 2451545.0D0 + ( DATE1 - 2000.D0 ) * 365.25D0
      ELSE
           TJD1 = DATE1
      END IF
      IF ( DATE2 .LT. 10000.D0 ) THEN
           TJD2 = 2451545.0D0 + ( DATE2 - 2000.D0 ) * 365.25D0
      ELSE
           TJD2 = DATE2
      END IF
 
* --- CONVERT INPUT ANGULAR COMPONENTS TO VECTORS ---------------------
 
*     IF PARALLAX IS UNKNOWN, UNDETERMINED, OR ZERO, SET IT TO 1E-6
*     MILLIARCSECOND, CORRESPONDING TO A DISTANCE OF 1 GIGAPARSEC
      PARALX = PARX1
      IF ( PARALX .LE. 0.D0 ) PARALX = 1.D-6
 
*     CONVERT RIGHT ASCENSION, DECLINATION, AND PARALLAX TO POSITION
*     VECTOR IN EQUATORIAL SYSTEM WITH UNITS OF AU
      DIST = 1.D0 / DSIN ( PARALX * 1.D-3 / SECCON )
      R = RA1 * 54000.D0 / SECCON
      D = DEC1 * 3600.D0 / SECCON
      CRA = DCOS ( R )
      SRA = DSIN ( R )
      CDC = DCOS ( D )
      SDC = DSIN ( D )
      POS1(1) = DIST * CDC * CRA
      POS1(2) = DIST * CDC * SRA
      POS1(3) = DIST * SDC
 
*     COMPUTE DOPPLER FACTOR, WHICH ACCOUNTS FOR CHANGE IN
*     LIGHT TRAVEL TIME TO STAR
      K = 1.D0 / ( 1.D0 - RV1 / C )
 
*     CONVERT PROPER MOTION AND RADIAL VELOCITY TO ORTHOGONAL
*     COMPONENTS OF MOTION, IN SPHERICAL POLAR SYSTEM AT STAR'S
*     ORIGINAL POSITION, WITH UNITS OF AU/DAY
      PMR = PMRA1  / ( PARALX * 365.25D0 ) * K
      PMD = PMDEC1 / ( PARALX * 365.25D0 ) * K
      RVL = RV1 * 86400.D0 / AUKM          * K
 
*     TRANSFORM MOTION VECTOR TO EQUATORIAL SYSTEM
      VEL1(1) = - PMR * SRA - PMD * SDC * CRA + RVL * CDC * CRA
      VEL1(2) =   PMR * CRA - PMD * SDC * SRA + RVL * CDC * SRA
      VEL1(3) =               PMD * CDC       + RVL * SDC
 
* --- UPDATE STAR'S POSITION VECTOR FOR SPACE MOTION ------------------
*     (ONLY IF IT=1 OR IT=3)
 
      IF ( IT .EQ. 1 .OR. IT .EQ. 3 ) THEN
          DO 22 J=1,3
              POS2(J) = POS1(J) + VEL1(J) * ( TJD2 - TJD1 )
              VEL2(J) = VEL1(J)
   22     CONTINUE
      ELSE
          DO 24 J=1,3
              POS2(J) = POS1(J)
              VEL2(J) = VEL1(J)
   24     CONTINUE
      END IF
 
* --- PRECESS POSITION AND VELOCITY VECTORS ---------------------------
*     (ONLY IF IT=2 OR IT=3)
 
      IF ( IT .EQ. 2 .OR. IT .EQ. 3 ) THEN
          DO 32 J=1,3
              POS1(J) = POS2(J)
              VEL1(J) = VEL2(J)
   32     CONTINUE
          CALL PRECES ( TJD1, POS1, TJD2,    POS2 )
          CALL PRECES ( TJD1, VEL1, TJD2,    VEL2 )
      END IF
 
* --- ROTATE DYNAMICAL J2000.0 POSITION AND VELOCITY VECTORS TO ICRS --
*     (ONLY IF IT=4)
 
      IF ( IT .EQ. 4 ) THEN
          CALL FRAME ( POS1, -1,    POS2 )
          CALL FRAME ( VEL1, -1,    VEL2 )
      END IF
 
* --- ROTATE ICRS POSITION AND VELOCITY VECTORS TO DYNAMICAL J2000.0 --
*     (ONLY IF IT=5)
 
      IF ( IT .EQ. 5 ) THEN
          CALL FRAME ( POS1, 1,    POS2 )
          CALL FRAME ( VEL1, 1,    VEL2 )
      END IF
 
* --- CONVERT VECTORS BACK TO ANGULAR COMPONENTS FOR OUTPUT -----------
 
*     FROM UPDATED POSITION VECTOR, OBTAIN STAR'S NEW POSITION
*     EXPRESSED AS ANGULAR QUANTITIES
      XYPROJ = DSQRT ( POS2(1)**2 + POS2(2)**2 )
      R = 0.D0
      IF ( XYPROJ .GT. 0.D0 ) R = DATAN2 ( POS2(2), POS2(1) )
      RA2 = R * SECCON / 54000.D0
      IF ( RA2 .LT.  0.D0 ) RA2 = RA2 + 24.D0
      IF ( RA2 .GE. 24.D0 ) RA2 = RA2 - 24.D0
      D = DATAN2 ( POS2(3), XYPROJ  )
      DEC2 = D * SECCON / 3600.D0
      DIST = DSQRT ( POS2(1)**2 + POS2(2)**2 + POS2(3)**2 )
      PARALX = DASIN ( 1.D0 / DIST ) * SECCON * 1.D3
      PARX2 = PARALX
 
*     TRANSFORM MOTION VECTOR BACK TO SPHERICAL POLAR SYSTEM AT STAR'S
*     NEW POSITION
      CRA = DCOS ( R )
      SRA = DSIN ( R )
      CDC = DCOS ( D )
      SDC = DSIN ( D )
      PMR = - VEL2(1) * SRA       + VEL2(2) * CRA
      PMD = - VEL2(1) * CRA * SDC - VEL2(2) * SRA * SDC + VEL2(3) * CDC
      RVL =   VEL2(1) * CRA * CDC + VEL2(2) * SRA * CDC + VEL2(3) * SDC
 
*     CONVERT COMPONENTS OF MOTION FROM AU/DAY TO NORMAL
*     CATALOG UNITS
      PMRA2  = PMR * PARALX * 365.25D0   / K
      PMDEC2 = PMD * PARALX * 365.25D0   / K
      RV2    = RVL * ( AUKM / 86400.D0 ) / K
 
*     TAKE CARE OF ZERO-PARALLAX CASE
      IF ( PARX2 .LE. 1.01D-6 ) THEN
          PARX2 = 0.D0
          RV2 = RV1
      END IF
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )
 
      RETURN
 
      END
 
 
 
      SUBROUTINE ZDAZ ( UJD, XP, YP, GLON, GLAT, HT, RA, DEC, IREFR,
     .                  ZD, AZ, RAR, DECR )
*
*     THIS SUBROUTINE TRANSFORMS TOPOCENTRIC RIGHT ASCENSION AND
*     DECLINATION TO ZENITH DISTANCE AND AZIMUTH.  THIS ROUTINE USES
*     A METHOD THAT PROPERLY ACCOUNTS FOR POLAR MOTION, WHICH IS
*     SIGNIFICANT AT THE SUB-ARCSECOND LEVEL.  THIS SUBROUTINE
*     CAN ALSO ADJUST COORDINATES FOR ATMOSPHERIC REFRACTION.
*
*          UJD    = UT1 JULIAN DATE (IN)
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
*                   IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF DATE (IN)
*          DEC    = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
*                   IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE (IN)
*          IREFR  = ATMOSPHERIC REFRACTION OPTION (IN)
*                   SET IREFR=0 FOR NO REFRACTION
*                   SET IREFR=1 TO INCLUDE REFRACTION
*          ZD     = TOPOCENTRIC ZENITH DISTANCE IN DEGREES,
*                   AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*          AZ     = TOPOCENTRIC AZIMUTH (MEASURED EAST FROM NORTH)
*                   IN DEGREES (OUT)
*          RAR    = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
*                   IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF DATE, AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*          DECR   = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
*                   IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE,
*                   AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*
*     NOTE 1:  XP AND YP CAN BE SET TO ZERO IF SUB-ARCSECOND ACCURACY IS
*     NOT NEEDED.  HT IS USED ONLY FOR REFRACTION, IF IREFR=1.  RA AND
*     DEC CAN BE OBTAINED FROM TPSTAR, TPPLAN, OR PLACE.
*
*     NOTE 2:  THE DIRECTONS ZD=0 (ZENITH) AND AZ=0 (NORTH) ARE
*     HERE CONSIDERED FIXED IN THE TERRESTRIAL SYSTEM.  SPECIFICALLY,
*     THE ZENITH IS ALONG THE GEODETIC NORMAL, AND NORTH IS TOWARD
*     THE ITRS POLE.
*
*     NOTE 3:  IF IREFR=0, THEN RAR=RA AND DECR=DEC.
*
*     NOTE 4: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION UJD,XP,YP,GLON,GLAT,HT,RA,DEC,ZD,AZ,RAR,DECR,
     .     PI,DEGRAD,RADDEG,
     .     SINLAT,COSLAT,SINLON,COSLON,SINDC,COSDC,SINRA,COSRA,
     .     UZE,UNE,UWE,UZ,UN,UW,P,PR,PZ,PN,PW,PROJ,
     .     ZD0,ZD1,REFR,SINZD,COSZD,SINZD0,COSZD0,
     .     DSIN,DCOS,DSQRT,DATAN2
      DIMENSION UZE(3), UNE(3), UWE(3), UZ(3), UN(3), UW(3),
     .     P(3), PR(3)
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           )
      PARAMETER ( RADDEG = 180.D0 / PI           )
 
      RAR    = RA
      DECR   = DEC
      SINLAT = DSIN ( GLAT * DEGRAD )
      COSLAT = DCOS ( GLAT * DEGRAD )
      SINLON = DSIN ( GLON * DEGRAD )
      COSLON = DCOS ( GLON * DEGRAD )
      SINDC  = DSIN ( DEC * DEGRAD )
      COSDC  = DCOS ( DEC * DEGRAD )
      SINRA  = DSIN ( RA * 15.D0 * DEGRAD )
      COSRA  = DCOS ( RA * 15.D0 * DEGRAD )
 
* --- SET UP ORTHONORMAL BASIS VECTORS IN LOCAL EARTH-FIXED SYSTEM ----
 
*     DEFINE VECTOR TOWARD LOCAL ZENITH IN EARTH-FIXED SYSTEM (Z AXIS)
      UZE(1) =  COSLAT * COSLON
      UZE(2) =  COSLAT * SINLON
      UZE(3) =  SINLAT
 
*     DEFINE VECTOR TOWARD LOCAL NORTH IN EARTH-FIXED SYSTEM (X AXIS)
      UNE(1) = -SINLAT * COSLON
      UNE(2) = -SINLAT * SINLON
      UNE(3) =  COSLAT
 
*     DEFINE VECTOR TOWARD LOCAL WEST IN EARTH-FIXED SYSTEM (Y AXIS)
      UWE(1) =  SINLON
      UWE(2) = -COSLON
      UWE(3) =  0.D0
 
* --- OBTAIN VECTORS IN CELESTIAL SYSTEM ------------------------------
 
*     ROTATE EARTH-FIXED ORTHONORMAL BASIS VECTORS TO CELESTIAL SYSTEM
*     (WRT EQUATOR AND EQUINOX OF DATE)
      CALL EQINOX
      CALL TERCEL ( -1.D0, UJD, XP, YP, UZE,   UZ )
      CALL TERCEL ( -1.D0, UJD, XP, YP, UNE,   UN )
      CALL TERCEL ( -1.D0, UJD, XP, YP, UWE,   UW )
      CALL RESUME
 
*     DEFINE UNIT VECTOR P TOWARD OBJECT IN CELESTIAL SYSTEM
*     (WRT EQUATOR AND EQUINOX OF DATE)
      P(1) = COSDC * COSRA
      P(2) = COSDC * SINRA
      P(3) = SINDC
 
* --- COMPUTE COORDINATES OF OBJECT WRT ORTHONORMAL BASIS -------------
 
*     COMPUTE COMPONENTS OF P -- PROJECTIONS OF P ONTO ROTATED
*     EARTH-FIXED BASIS VECTORS
      PZ = P(1) * UZ(1) + P(2) * UZ(2) + P(3) * UZ(3)
      PN = P(1) * UN(1) + P(2) * UN(2) + P(3) * UN(3)
      PW = P(1) * UW(1) + P(2) * UW(2) + P(3) * UW(3)
 
*     COMPUTE AZIMUTH AND ZENITH DISTANCE
      PROJ = DSQRT ( PN**2 + PW**2 )
      AZ = 0.D0
      IF ( PROJ .GT. 0.D0 ) AZ = -DATAN2 ( PW, PN ) * RADDEG
      IF ( AZ .LT.   0.D0 ) AZ = AZ + 360.D0
      IF ( AZ .GE. 360.D0 ) AZ = AZ - 360.D0
      ZD = DATAN2 ( PROJ, PZ ) * RADDEG
 
* --- APPLY ATMOSPHERIC REFRACTION IF REQUESTED -----------------------
 
      IF ( IREFR .EQ. 1 ) THEN
 
*         GET REFRACTION IN ZENITH DISTANCE
*         ITERATIVE PROCESS REQUIRED BECAUSE REFRACTION ALGORITHMS ARE
*         ALWAYS A FUNCTION OF OBSERVED (NOT COMPUTED) ZENITH DISTANCE
          ZD0 = ZD
  40      ZD1 = ZD
          CALL REFRAC ( HT, ZD,   REFR )
          ZD = ZD0 - REFR
*         REQUIRE CONVERGENCE TO 0.1 ARCSEC (ACTUAL ACCURACY LESS)
          IF ( DABS ( ZD - ZD1 ) .GT. 3.D-5 ) GO TO 40
 
*         APPLY REFRACTION TO CELESTIAL COORDINATES OF OBJECT
          IF ( REFR .GT. 0.D0 .AND. ZD .GT. 3.D-4 ) THEN
 
*             SHIFT POSITION VECTOR OF OBJECT IN CELESTIAL SYSTEM
*             TO ACCOUNT FOR FOR REFRACTION (SEE USNO/AA TECHNICAL
*             NOTE 1998-09)
              SINZD  = DSIN ( ZD * DEGRAD )
              COSZD  = DCOS ( ZD * DEGRAD )
              SINZD0 = DSIN ( ZD0 * DEGRAD )
              COSZD0 = DCOS ( ZD0 * DEGRAD )
*             COMPUTE REFRACTED POSITION VECTOR
              DO J = 1, 3
  50          PR(J) = ( ( P(J) - COSZD0 * UZ(J) ) / SINZD0 ) * SINZD
     .                 +                  UZ(J)              * COSZD
                enddo
*             COMPUTE REFRACTED RIGHT ASCENSION AND DECLINATION
              PROJ = DSQRT ( PR(1)**2 + PR(2)**2 )
              RAR = 0.D0
              IF ( PROJ .GT. 0.D0 ) RAR = DATAN2 ( PR(2), PR(1) )
     .                                    * RADDEG / 15.D0
              IF ( RAR .LT.  0.D0 ) RAR = RAR + 24.D0
              IF ( RAR .GE. 24.D0 ) RAR = RAR - 24.D0
              DECR = DATAN2 ( PR(3), PROJ ) * RADDEG
 
          END IF
 
      END IF
 
* ---------------------------------------------------------------------
 
      RETURN
 
      END
 
 
 
      SUBROUTINE GCRSEQ ( TJD, ICOORD, RAG, DECG,   RA, DEC )
*
*     THIS SUBROUTINE CONVERTS GCRS RIGHT ASCENSION AND DECLINATION
*     TO COORDINATES WITH RESPECT TO THE EQUATOR OF DATE (MEAN OR TRUE).
*     FOR COORDINATES WITH RESPECT TO THE TRUE EQUATOR OF DATE, THE
*     ORIGIN OF RIGHT ASCENSION CAN BE EITHER THE TRUE EQUINOX OR THE
*     CELESTIAL INTERMEDIATE ORIGIN (CIO).
*
*          TJD    = TT JULIAN DATE OF EQUATOR TO BE USED FOR
*                   OUTPUT COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION FOR OUTPUT
*                   COORDINATES (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=2 FOR TRUE EQUATOR AND CIO OF DATE
*          RAG    = GCRS RIGHT ASCENSION IN HOURS (IN)
*          DECG   = GCRS DECLINATION IN DEGREES (IN)
*          RA     = RIGHT ASCENSION IN HOURS, REFERRED TO SPECIFIED
*                   EQUATOR AND RIGHT ASCENSION ORIGIN OF DATE (OUT)
*          DEC    = DECLINATION IN DEGREES, REFERRED TO SPECIFIED
*                   EQUATOR OF DATE (OUT)
*
*
      DOUBLE PRECISION TJD,RAG,DECG,RA,DEC,PI,RADCON,T0,T1,T,SECDIF,R,D,
     .     POS1,POS2,POS3,POS4,RCIO,X,Y,Z,DSIN,DCOS
      DIMENSION POS1(3), POS2(3), POS3(3), POS4(3), X(3), Y(3), Z(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
 
*     T1 IS THE TDB JULIAN DATE
      CALL TIMES ( TJD, T,   SECDIF )
      T1 = TJD + SECDIF / 86400.D0
 
*     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
      R = RAG * 15.D0 * RADCON
      D = DECG * RADCON
      POS1(1) = DCOS ( D ) * DCOS ( R )
      POS1(2) = DCOS ( D ) * DSIN ( R )
      POS1(3) = DSIN ( D )
 
      IF ( ICOORD .LE. 1 ) THEN
 
*         TRANSFORM THE POSITION VECTOR FROM GCRS TO MEAN EQUATOR AND
*         EQUINOX OF DATE
          CALL FRAME  ( POS1, 1,   POS2 )
          CALL PRECES ( T0, POS2, T1,   POS3 )
*         IF REQUESTED, TRANSFORM FURTHER TO TRUE EQUATOR AND EQUINOX
*         OF DATE
          IF ( ICOORD .EQ. 1 ) THEN
              CALL NUTATE ( T1, POS3,   POS4)
          ELSE
              POS4(1) = POS3(1)
              POS4(2) = POS3(2)
              POS4(3) = POS3(3)
          END IF
 
      ELSE
 
*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( T1,   RCIO, KCIO )
          CALL CIOBAS ( T1, RCIO, KCIO,   X, Y, Z )
 
*         TRANSFORM POSITION VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
*         (WHICH HAS THE CIO AS ITS ORIGIN OF RIGHT ASCENSION)
          POS4(1) = X(1) * POS1(1) + X(2) * POS1(2) + X(3) * POS1(3)
          POS4(2) = Y(1) * POS1(1) + Y(2) * POS1(2) + Y(3) * POS1(3)
          POS4(3) = Z(1) * POS1(1) + Z(2) * POS1(2) + Z(3) * POS1(3)
 
      END IF
 
      CALL ANGLES ( POS4,   RA, DEC )
 
      RETURN
 
      END
 
 
 
      SUBROUTINE EQECL ( TJD, ICOORD, RA, DEC,   ELON, ELAT )
*
*     THIS SUBROUTINE CONVERTS RIGHT ASCENSION AND DECLINATION
*     TO ECLIPTIC LONGITUDE AND LATITUDE.
*
*          TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
*                   USED FOR COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   (ECLIPTIC IS ALWAYS THE MEAN PLANE)
*          RA     = RIGHT ASCENSION IN HOURS, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (IN)
*          DEC    = DECLINATION IN DEGREES, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (IN)
*          ELON   = ECLIPTIC LONGITUDE IN DEGREES, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (OUT)
*          ELAT   = ECLIPTIC LATITUDE IN DEGREES, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (OUT)
*
*     NOTE:  TO CONVERT ICRS RA AND DEC TO ECLIPTIC COORDINATES (MEAN
*     ECLIPTIC AND EQUINOX OF J2000.0), SET TJD = 0.D0 AND ICOORD = 0.
*     EXCEPT FOR THE INPUT TO THIS CASE, ALL COORDINATES ARE DYNAMICAL.
*
*
      DOUBLE PRECISION TJD,RA,DEC,ELON,ELAT,PI,RADCON,R,D,XYPROJ,E,
     .     POS1,POS2,DSIN,DCOS,DSQRT,DATAN2
      DIMENSION POS1(3), POS2(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
*     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
      R = RA * 15.D0 * RADCON
      D = DEC * RADCON
      POS1(1) = DCOS ( D ) * DCOS ( R )
      POS1(2) = DCOS ( D ) * DSIN ( R )
      POS1(3) = DSIN ( D )
 
*     CONVERT THE VECTOR FROM EQUATORIAL TO ECLIPTIC SYSTEM
      CALL EQEC ( TJD, ICOORD, POS1,   POS2 )
 
*     DECOMPOSE ECLIPTIC VECTOR INTO ECLIPTIC LONGITUDE AND LATITUDE
      XYPROJ = DSQRT ( POS2(1)**2 + POS2(2)**2 )
      E = 0.D0
      IF ( XYPROJ .GT. 0.D0 ) E = DATAN2 ( POS2(2), POS2(1) )
      ELON = E / RADCON
      IF ( ELON .LT. 0.D0 ) ELON = ELON + 360.D0
      E = DATAN2 ( POS2(3), XYPROJ )
      ELAT = E / RADCON
 
      RETURN
 
      END
 
 
 
      SUBROUTINE EQEC ( TJD, ICOORD, POS1,   POS2 )
*
*     THIS SUBROUTINE CONVERTS AN EQUATORIAL POSITION VECTOR TO
*     AN ECLIPTIC POSITION VECTOR.
*
*          TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
*                   USED FOR COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   (ECLIPTIC IS ALWAYS THE MEAN PLANE)
*          POS1   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (IN)
*          POS2   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (OUT)
*
*     NOTE:  TO CONVERT ICRS VECTORS TO ECLIPTIC VECTORS (MEAN ECLIPTIC
*     AND EQUINOX OF J2000.0 ONLY), SET TJD = 0.D0 AND ICOORD = 0.
*     EXCEPT FOR THE INPUT TO THIS CASE, ALL VECTORS ARE ASSUMED TO
*     BE WITH RESPECT TO A DYNAMICAL SYSTEM.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,POS0,PI,RADCON,T0,T1,T,SECDIF,
     .     TLAST,OB2000,OBLM,OBLT,OBL,X,DSIN,DCOS
      DIMENSION POS1(3), POS2(3), POS0(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   OB2000 / 0.D0 /
 
*     T1 IS THE TDB JULIAN DATE
      CALL TIMES ( TJD, T,   SECDIF )
      T1 = TJD + SECDIF / 86400.D0
 
      IF ( TJD .EQ. 0.D0 ) THEN
*         CASE WHERE INPUT VECTOR IS IN ICRS SYSTEM
          CALL FRAME ( POS1, 1,   POS0 )
*         GET MEAN OBLIQUITY AT J2000.0 IF NECESSARY
          IF ( OB2000 .EQ. 0.D0 ) CALL ETILT ( T0,  OB2000, X, X, X, X )
          OBL = OB2000 * RADCON
      ELSE
*         CASE WHERE INPUT VECTOR IS IN EQUATOR OF DATE SYSTEM
          POS0(1) = POS1(1)
          POS0(2) = POS1(2)
          POS0(3) = POS1(3)
*         GET MEAN AND TRUE OBLIQUITY
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( T1,   OBLM, OBLT, X, X, X )
              TLAST = TJD
          END IF
*         SELECT MEAN OR TRUE OBLIQUITY
          OBL = OBLM * RADCON
          IF ( ICOORD .EQ. 1 ) OBL = OBLT * RADCON
      END IF
 
*     ROTATE EQUATORIAL POSITION VECTOR TO ECLIPTIC SYSTEM
      POS2(1) =  POS0(1)
      POS2(2) =  POS0(2) * DCOS ( OBL ) + POS0(3) * DSIN ( OBL )
      POS2(3) = -POS0(2) * DSIN ( OBL ) + POS0(3) * DCOS ( OBL )
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )
 
      RETURN
 
      END
 
 
 
      SUBROUTINE ECEQ ( TJD, ICOORD, POS1,   POS2 )
*
*     THIS SUBROUTINE CONVERTS AN ECLIPTIC POSITION VECTOR TO
*     AN EQUATORIAL POSITION VECTOR.
*
*          TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
*                   USED FOR COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   (ECLIPTIC IS ALWAYS THE MEAN PLANE)
*          POS1   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (IN)
*          POS2   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (OUT)
*
*     NOTE:  TO CONVERT ECLIPTIC VECTORS (MEAN ECLIPTIC AND EQUINOX OF
*     OF J2000.0 ONLY) TO ICRS VECTORS, SET TJD = 0.D0 AND ICOORD = 0.
*     EXCEPT FOR THE OUTPUT FROM THIS CASE, ALL VECTORS ARE ASSUMED TO
*     BE WITH RESPECT TO A DYNAMICAL SYSTEM.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,POS0,PI,RADCON,T0,T1,T,SECDIF,
     .     TLAST,OB2000,OBLM,OBLT,OBL,X,DSIN,DCOS
      DIMENSION POS1(3), POS2(3), POS0(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   OB2000 / 0.D0 /
 
*     T1 IS THE TDB JULIAN DATE
      CALL TIMES ( TJD, T,   SECDIF )
      T1 = TJD + SECDIF / 86400.D0
 
      IF ( TJD .EQ. 0.D0 ) THEN
*         CASE WHERE OUTPUT VECTOR IS TO BE IN ICRS SYSTEM
*         GET MEAN OBLIQUITY AT J2000.0 IF NECESSARY
          IF ( OB2000 .EQ. 0.D0 ) CALL ETILT ( T0,  OB2000, X, X, X, X )
          OBL = OB2000 * RADCON
      ELSE
*         CASE WHERE OUTPUT VECTOR IS TO BE IN EQUATOR OF DATE SYSTEM
*         GET MEAN AND TRUE OBLIQUITY
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( T1,   OBLM, OBLT, X, X, X )
              TLAST = TJD
          END IF
*         SELECT MEAN OR TRUE OBLIQUITY
          OBL = OBLM * RADCON
          IF ( ICOORD .EQ. 1 ) OBL = OBLT * RADCON
      END IF
 
*     ROTATE ECLIPTIC POSITION VECTOR TO EQUATORIAL SYSTEM
      POS2(1) =  POS1(1)
      POS2(2) =  POS1(2) * DCOS ( OBL ) - POS1(3) * DSIN ( OBL )
      POS2(3) =  POS1(2) * DSIN ( OBL ) + POS1(3) * DCOS ( OBL )
 
      IF ( TJD .EQ. 0.D0 ) THEN
*         CASE WHERE OUTPUT VECTOR IS TO BE IN ICRS SYSTEM
          POS0(1) = POS2(1)
          POS0(2) = POS2(2)
          POS0(3) = POS2(3)
          CALL FRAME ( POS0, -1,   POS2 )
      END IF
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )
 
      RETURN
 
      END
 
 
 
      SUBROUTINE EQGAL ( RA, DEC,   GLON, GLAT )
*
*     THIS SUBROUTINE CONVERTS ICRS RIGHT ASCENSION AND DECLINATION
*     TO GALACTIC LONGITUDE AND LATITUDE.  IT USES THE TRANSFORMATION
*     GIVEN IN THE HIPPARCOS AND TYCHO CATALOGUES, VOL. 1,
*     SECTION 1.5.3.
*
*          RA     = ICRS RIGHT ASCENSION IN HOURS (IN)
*          DEC    = ICRS DECLINATION IN DEGREES (IN)
*          GLON   = GALACTIC LONGITUDE IN DEGREES (OUT)
*          GLAT   = GALACTIC LATITUDE IN DEGREES (OUT)
*
*
      DOUBLE PRECISION RA,DEC,GLON,GLAT,PI,RADCON,AG,R,D,XYPROJ,G,
     .     POS1,POS2,DSIN,DCOS,DSQRT,DATAN2
      DIMENSION POS1(3), POS2(3), AG(3,3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
*     ROTATION MATRIX A_G FROM HIPPARCOS DOCUMENTATION EQ. 1.5.11
      DATA AG /
     .     -0.0548755604D0, +0.4941094279D0, -0.8676661490D0,
     .     -0.8734370902D0, -0.4448296300D0, -0.1980763734D0,
     .     -0.4838350155D0, +0.7469822445D0, +0.4559837762D0 /
 
*     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
      R = RA * 15.D0 * RADCON
      D = DEC * RADCON
      POS1(1) = DCOS ( D ) * DCOS ( R )
      POS1(2) = DCOS ( D ) * DSIN ( R )
      POS1(3) = DSIN ( D )
 
*     ROTATE POSITION VECTOR TO GALACTIC SYSTEM, USING HIPPARCOS
*     DOCUMENTATION EQ. 1.5.13
      POS2(1) = AG(1,1)*POS1(1) + AG(1,2)*POS1(2) + AG(1,3)*POS1(3)
      POS2(2) = AG(2,1)*POS1(1) + AG(2,2)*POS1(2) + AG(2,3)*POS1(3)
      POS2(3) = AG(3,1)*POS1(1) + AG(3,2)*POS1(2) + AG(3,3)*POS1(3)
 
*     DECOMPOSE GALACTIC VECTOR INTO LONGITUDE AND LATITUDE
      XYPROJ = DSQRT ( POS2(1)**2 + POS2(2)**2 )
      G = 0.D0
      IF ( XYPROJ .GT. 0.D0 ) G = DATAN2 ( POS2(2), POS2(1) )
      GLON = G / RADCON
      IF ( GLON .LT. 0.D0 ) GLON = GLON + 360.D0
      G = DATAN2 ( POS2(3), XYPROJ )
      GLAT = G / RADCON
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )
 
      RETURN
 
      END
 
 
 
      SUBROUTINE VECTRS (RA,DEC,PMRA,PMDEC,PARLLX,RV,POS,VEL)
*
*     THIS SUBROUTINE CONVERTS ANGULAR QUANTITIES RELATED TO A STAR'S
*     POSITION AND MOTION TO VECTORS.
*
*          RA     = RIGHT ASCENSION IN HOURS (IN)
*          DEC    = DECLINATION IN DEGREES (IN)
*          PMRA   = PROPER MOTION IN RA IN MILLIARCSECONDS PER YEAR
*                   (IN)
*          PMDEC  = PROPER MOTION IN DEC IN MILLIARCSECONDS PER YEAR
*                   (IN)
*          PARLLX = PARALLAX IN MILLIARCSECONDS (IN)
*          RV     = RADIAL VELOCITY IN KILOMETERS/SECOND (IN)
*          POS    = POSITION VECTOR, EQUATORIAL RECTANGULAR COORDINATES,
*                   WITH RESPECT TO SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU (OUT)
*          VEL    = VELOCITY VECTOR, EQUATORIAL RECTANGULAR COORDINATES,
*                   WITH RESPECT TO SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU/DAY (OUT)
*
*
      DOUBLE PRECISION RA,DEC,PMRA,PMDEC,PARLLX,RV,POS,VEL,
     .     PI,SECCON,AUKM,C,PARALX,DIST,R,D,CRA,SRA,CDC,SDC,K,
     .     PMR,PMD,RVL,DCOS,DSIN
      DIMENSION POS(3), VEL(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET LENGTH OF AU IN KILOMETERS
          CALL ASTCON ( 'AU', 1.D-3,   AUKM )
*         GET C, THE SPEED OF LIGHT IN KILOMETERS/SECOND
          CALL ASTCON ( 'C', 1.D-3,   C )
      END IF
 
*     IF PARALLAX IS UNKNOWN, UNDETERMINED, OR ZERO, SET IT TO 1E-6
*     MILLIARCSECOND, CORRESPONDING TO A DISTANCE OF 1 GIGAPARSEC
      PARALX = PARLLX
      IF ( PARALX .LE. 0.D0 ) PARALX = 1.D-6
 
*     CONVERT RIGHT ASCENSION, DECLINATION, AND PARALLAX TO POSITION
*     VECTOR IN EQUATORIAL SYSTEM WITH UNITS OF AU
      DIST = 1.D0 / DSIN ( PARALX * 1.D-3 / SECCON )
      R = RA * 54000.D0 / SECCON
      D = DEC * 3600.D0 / SECCON
      CRA = DCOS ( R )
      SRA = DSIN ( R )
      CDC = DCOS ( D )
      SDC = DSIN ( D )
      POS(1) = DIST * CDC * CRA
      POS(2) = DIST * CDC * SRA
      POS(3) = DIST * SDC
 
*     COMPUTE DOPPLER FACTOR, WHICH ACCOUNTS FOR CHANGE IN
*     LIGHT TRAVEL TIME TO STAR
      K = 1.D0 / ( 1.D0 - RV / C )
 
*     CONVERT PROPER MOTION AND RADIAL VELOCITY TO ORTHOGONAL COMPONENTS
*     OF MOTION WITH UNITS OF AU/DAY
      PMR = PMRA  / ( PARALX * 365.25D0 ) * K
      PMD = PMDEC / ( PARALX * 365.25D0 ) * K
      RVL = RV * 86400.D0 / AUKM          * K
 
*     TRANSFORM MOTION VECTOR TO EQUATORIAL SYSTEM
      VEL(1) = - PMR * SRA - PMD * SDC * CRA + RVL * CDC * CRA
      VEL(2) =   PMR * CRA - PMD * SDC * SRA + RVL * CDC * SRA
      VEL(3) =               PMD * CDC       + RVL * SDC
 
      RETURN
 
      END
 
 
 
      SUBROUTINE ANGLES (POS,RA,DEC)
*
*     THIS SUBROUTINE CONVERTS A VECTOR TO ANGULAR QUANTITIES.
*
*          POS = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                COORDINATES (IN)
*          RA  = RIGHT ASCENSION IN HOURS (OUT)
*          DEC = DECLINATION IN DEGREES (OUT)
*
*
      DOUBLE PRECISION POS,RA,DEC,PI,SECCON,XYPROJ,R,D,DSQRT,DATAN2
      DIMENSION POS(3)
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
      XYPROJ = DSQRT(POS(1)**2 + POS(2)**2)
      R = 0.D0
      IF (XYPROJ.GT.0.D0) R = DATAN2(POS(2),POS(1))
      RA = R * SECCON / 54000.D0
      IF (RA.LT. 0.D0) RA = RA + 24.D0
      IF (RA.GE.24.D0) RA = RA - 24.D0
      D = DATAN2(POS(3),XYPROJ)
      DEC = D * SECCON / 3600.D0
 
*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC (POS)
 
      RETURN
 
      END
 
 
 
      SUBROUTINE PROPMO (TJD1,POS1,VEL1,TJD2,POS2)
*
*     THIS SUBROUTINE APPLIES PROPER MOTION, INCLUDING FORESHORTENING
*     EFFECTS, TO A STAR'S POSITION.
*
*          TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
*          POS1 = POSITION VECTOR OF STAR AT FIRST EPOCH (IN)
*          VEL1 = VELOCITY VECTOR OF STAR AT FIRST EPOCH (IN)
*          TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
*          POS2 = POSITION VECTOR OF STAR AT SECOND EPOCH (OUT)
*
*
      DOUBLE PRECISION TJD1,POS1,VEL1,TJD2,POS2
      DIMENSION POS1(3), VEL1(3), POS2(3)
 
      DO J=1,3
   20 POS2(J) = POS1(J) + VEL1(J) * (TJD2 - TJD1)
         enddo
      RETURN
 
      END
 
 
 
      SUBROUTINE GEOCEN (POS1,PE,POS2,TLIGHT)
*
*     THIS SUBROUTINE MOVES THE ORIGIN OF COORDINATES FROM THE
*     BARYCENTER OF THE SOLAR SYSTEM TO THE OBSERVER (OR THE
*     GEOCENTER).  I.E., THIS SUBROUTINE ACCOUNTS FOR PARALLAX
*     (ANNUAL+GEOCENTRIC OR JUST ANNUAL).
*
*          POS1   = POSITION VECTOR OF STAR OR PLANET, WITH RESPECT TO
*                   ORIGIN AT SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU (IN)
*          PE     = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU (IN)
*          POS2   = POSITION VECTOR OF STAR OR PLANET, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (OUT)
*          TLIGHT = LIGHT-TIME FROM STAR OR PLANET TO OBSERVER (OR THE
*                   GEOCENTER) IN DAYS (OUT)
*
*     NOTE: STAR AND PLANET ARE USED GENERICALLY FOR BODIES OUTSIDE AND
*           INSIDE THE SOLAR SYSTEM, RESPECTIVELY.
*
      DOUBLE PRECISION POS1,PE,POS2,TLIGHT,C,DSQRT
      DIMENSION POS1(3), PE(3), POS2(3)
      SAVE
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ('C(AU/DAY)',1.D0,C)
      END IF
 
      DO J=1,3
   20 POS2(J) = POS1(J) - PE(J)
         enddo
      TLIGHT = DSQRT(POS2(1)**2 + POS2(2)**2 + POS2(3)**2) / C
 
      RETURN
 
      END
 
 
 
      SUBROUTINE GEOPOS (TJD,LOCATN,OBSERV,POS,VEL)
*
*     THIS SUBROUTINE COMPUTES THE GEOCENTRIC POSITION AND VELOCITY
*     OF AN OBSERVER ON THE SURFACE OF THE EARTH OR ON A NEAR-EARTH
*     SPACECRAFT.  THE FINAL VECTORS ARE EXPRESSED IN THE GCRS.
*
*          TJD    = TT JULIAN DATE (IN)
*          LOCATN = INTEGER CODE SPECIFYING LOCATION OF OBSERVER (IN)
*                   SET LOCATN=0 FOR OBSERVER AT GEOCENTER
*                   SET LOCATN=1 FOR OBSERVER ON SURFACE OF EARTH
*                   SET LOCATN=2 FOR OBSERVER ON NEAR-EARTH SPACECRAFT
*          OBSERV = ARRAY OF DATA SPECIFYING LOCATION OF OBSERVER (IN)
*                   FOR LOCATN=0, THIS ARRAY NOT USED
*                   FOR LOCATN=1,
*                   OBSERV(1) = GEODETIC LONGITUDE (WGS-84) OF OBSERVER
*                               (EAST +) IN DEGREES (IN)
*                   OBSERV(2) = GEODETIC LATITUDE (WGS-84) OF OBSERVER
*                               (NORTH +) IN DEGREES (IN)
*                   OBSERV(3) = HEIGHT OF OBSERVER ABOVE ELLIPSOID
*                               IN METERS (IN)
*                   OBSERV(4) = VALUE OF DELTA-T IN SECONDS (IN)
*                               (DELTA-T=TT-UT1)
*                   OBSERV(5) = (NOT USED, RESERVED FOR FUTURE USE)
*                   OBSERV(6) = (NOT USED, RESERVED FOR FUTURE USE)
*                   FOR LOCATN=2,
*                   OBSERV(1) = GEOCENTRIC X IN KILOMETERS
*                   OBSERV(2) = GEOCENTRIC Y IN KILOMETERS
*                   OBSERV(3) = GEOCENTRIC Z IN KILOMETERS
*                   OBSERV(4) = GEOCENTRIC X-DOT IN KILOMETERS/SECOND
*                   OBSERV(5) = GEOCENTRIC Y-DOT IN KILOMETERS/SECOND
*                   OBSERV(6) = GEOCENTRIC Z-DOT IN KILOMETERS/SECOND
*                   WITH RESPECT TO TRUE EQUATOR AND EQUINOX OF DATE
*          POS    = POSITION VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
*                   AT GEOCENTER, REFERRED TO GCRS AXES, COMPONENTS
*                   IN AU (OUT)
*          VEL    = VELOCITY VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
*                   AT GEOCENTER, REFERRED TO GCRS AXES, COMPONENTS
*                   IN AU/DAY (OUT)
*
*     NOTE 1: IF LOCATN=1 AND OBSERV(4)=0.D0, THE VALUE OF DELTA-T WILL
*     BE OBTAINED FROM GETDT, WHICH PROVIDES THE LAST VALUE OF DELTA-T
*     DEFINED BY USER VIA CALL TO SETDT.
*
*     NOTE 2: THIS SUBROUTINE CALLS SUBROUTINE TERRA FOR AN OBSERVER
*     ON THE SURFACE OF THE EARTH.  TERRA NEGLECTS POLAR MOTION, AN
*     APPROXIMATION WHICH MAY YIELD UP TO 15 METERS ERROR IN POSITION
*     AND SEVERAL MILLIMETERS/SEC ERROR IN VELOCITY.
*
*
      DOUBLE PRECISION TJD,OBSERV,POS,VEL,T0,TLAST,
     .     AU,DELTAT,TTJD,TDBJD,UT1JD,ST,GST,GMST,GAST,EQEQ,X,
     .     POS1,VEL1,POS2,VEL2,POS3,VEL3
 
      DIMENSION OBSERV(6), POS(3), VEL(3),
     .     POS1(3), VEL1(3), POS2(3), VEL2(3), POS3(3), VEL3(3)
 
      SAVE
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   GST / -99.D0 /,   NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET AU, THE LENGTH OF THE ASTRONOMICAL UNIT IN KILOMETERS
          CALL ASTCON ( 'AU', 1.D-3,   AU )
      END IF
 
      IF ( LOCATN .EQ. 0 ) THEN
          POS(1) = 0.D0
          POS(2) = 0.D0
          POS(3) = 0.D0
          VEL(1) = 0.D0
          VEL(2) = 0.D0
          VEL(3) = 0.D0
          GO TO 70
      END IF
 
      TTJD  = TJD
*     TDB IS APPROXIMATED BY TT
      TDBJD = TJD
 
*     GET GEOCENTRIC POSITION AND VELOCITY VECTORS OF OBSERVER WRT
*     EQUATOR AND EQUINOX OF DATE
 
      IF ( LOCATN .EQ. 1 ) THEN
 
*         OBSERVER ON SURFACE OF EARTH
 
*         TEMPORARY CODE TO USE SIDEREAL TIME PREVIOUSLY PROVIDED
          IF ( GST .NE. -99.D0 ) THEN
              GAST = GST
              GST = -99.D0
              GO TO 20
          END IF
*         END OF TEMPROARY CODE
 
*         GET DELTA-T VALUE
          IF ( OBSERV(4) .NE. 0.D0 ) THEN
              DELTAT = OBSERV(4) / 86400.D0
          ELSE
              CALL GETDT ( DELTAT )
          END IF
 
*         USING DELTA-T VALUE, COMPUTE UT1 AND SIDEREAL TIME
          IF ( TTJD .EQ. 0.D0 ) THEN
              UT1JD = TDBJD - DELTAT
          ELSE
              UT1JD = TTJD - DELTAT
          END IF
          IF ( DABS ( UT1JD - TLAST ) .GT. 1.D-8 ) THEN
              CALL EQINOX
              CALL SIDTIM ( UT1JD, 0.D0, 0,   GMST )
              CALL ETILT ( TDBJD,   X, X, EQEQ, X, X )
              CALL RESUME
              TLAST = UT1JD
          END IF
          GAST = GMST + EQEQ / 3600.D0
 
*         SUBROUTINE TERRA DOES THE HARD WORK, GIVEN SIDEREAL TIME
  20      CALL TERRA ( OBSERV(1), OBSERV(2), OBSERV(3), GAST,
     .                 POS1, VEL1 )
 
 
      ELSE IF ( LOCATN .EQ. 2 ) THEN
 
*         OBSERVER ON NEAR-EARTH SPACECRAFT
 
*         CONVERT UNITS TO AU AND AU/DAY
          DO 30 J = 1, 3
              POS1(J) = OBSERV(J)   / AU
              VEL1(J) = OBSERV(J+3) / AU * 86400.D0
  30      CONTINUE
 
      END IF
 
*     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
      CALL NUTATE ( -TDBJD, POS1,   POS2 )
      CALL PRECES ( TDBJD, POS2, T0,   POS3 )
      CALL FRAME ( POS3, -1,   POS )
 
*     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
      CALL NUTATE ( -TDBJD, VEL1,   VEL2 )
      CALL PRECES ( TDBJD, VEL2, T0,   VEL3 )
      CALL FRAME ( VEL3, -1,   VEL )
 
  70  RETURN
 
 
*     TEMPORARY CODE FOR COMPATIBILITY WITH OLD ROUTINES
      ENTRY PLACST ( ST )
      GST = ST
      RETURN
 
      END
 
 
 
      SUBROUTINE LITTIM (TJD,IDBODY,POSE,TLITE,POS,TLIGHT)
*
*     THIS SUBROUTINE COMPUTES THE POSITION OF A SOLAR SYSTEM BODY,
*     AS ANTEDATED FOR LIGHT-TIME.
*
*          TJD    = TDB JULIAN DATE OF OBSERVATION (IN)
*          IDBODY = ID NUMBER OF BODY, USED IN CALLS TO SOLSYS (IN)
*          POSE   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   REFERRED TO ICRS AXES, COMPONENTS IN AU (IN)
*          TLITE  = FIRST APPROXIMATION TO LIGHT-TIME, IN DAYS (IN)
*                   (CAN BE SET TO 0.D0 IF UNKNOWN)
*          POS    = POSITION VECTOR OF BODY, WITH RESPECT TO ORIGIN AT
*                   OBSERVER (OR THE GEOCENTER), REFERRED TO ICRS AXES,
*                   COMPONENTS IN AU (OUT)
*          TLIGHT = FINAL LIGHT-TIME, IN DAYS (OUT)
*
*
      DOUBLE PRECISION TJD,POSE,TLITE,POS,TLIGHT,T0,T1,T2,T3,TOL,
     .     POS1,VEL1,DINT,DABS
      LOGICAL SPLIT
 
      DIMENSION POSE(3), POS(3), POS1(3), VEL1(3)
 
      SAVE NTIMES, SPLIT
 
      DATA NTIMES / 0 /,   SPLIT / .FALSE. /
 
   3  FORMAT ( ' LITTIM: PROBLEM WITH BODY NUMBER ', I3, ' AT JD ',
     .     F10.1 )
 
      NTIMES = NTIMES + 1
 
*     ON FIRST CALL, CHECK WHETHER SOLSYS SUPPORTS SPLIT JULIAN DATES
      IF ( NTIMES .EQ. 1 ) SPLIT = IDSS ('JD') .EQ. 2
 
*     SET LIGHT-TIME CONVERGENCE TOLERANCE
      TOL = 1.D-9
      IF ( SPLIT .AND. TLITE .LT. 0.01D0 ) TOL = 1.D-12
 
*     IF SOLSYS SUPPORTS SPLIT JULIAN DATES, SPLIT THE JULIAN DATE
*     INTO WHOLE DAYS + FRACTION OF DAY
      T0 = 0.D0
      IF ( SPLIT ) T0 = DINT ( TJD )
      T1 = TJD - T0
      T2 = T1 - TLITE
      IF ( SPLIT ) CALL SOLSYS ( T0, IDBODY, 0,   POS1, VEL1, IERR)
      ITER = 0
 
*     ITERATE TO OBTAIN CORRECT LIGHT-TIME (USUALLY CONVERGES RAPIDLY)
  40  CALL SOLSYS ( T2, IDBODY, 0,   POS1, VEL1, IERR )
      CALL GEOCEN ( POS1, POSE,   POS, TLIGHT )
      IF ( IERR .NE. 0 ) THEN
          WRITE ( *, 3 ) IDBODY, T0 + T2
          GO TO 70
      END IF
      T3 = T1 - TLIGHT
      IF ( DABS ( T3 - T2 ) .GT. TOL ) THEN
          ITER = ITER + 1
          IF ( ITER .GT. 10 ) THEN
              WRITE ( *, 3 ) IDBODY, T0 + T3
              GO TO 70
          END IF
          T2 = T3
          GO TO 40
      END IF
 
  70  RETURN
 
      END
 
 
 
      SUBROUTINE DLIGHT (POS1,PE,DIFLT)
*
*     THIS SUBROUTINE RETURNS THE DIFFERENCE IN LIGHT-TIME, FOR A STAR,
*     BETWEEN THE BARYCENTER OF THE SOLAR SYSTEM AND THE OBSERVER (OR
*     THE GEOCENTER).
*
*          POS1   = POSITION VECTOR OF STAR, WITH RESPECT TO ORIGIN AT
*                   SOLAR SYSTEM BARYCENTER (IN)
*          PE     = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU (IN)
*          DIFLT  = DIFFERENCE IN LIGHT TIME, IN THE SENSE STAR TO
*                   BARYCENTER MINUS STAR TO EARTH, IN DAYS (OUT)
*
*     -OR-
*
*     THIS SUBROUTINE RETURNS THE LIGHT-TIME FROM THE OBSERVER (OR THE
*     GEOCENTER) TO A POINT ON A LIGHT RAY THAT IS CLOSEST TO A
*     SPECIFIC SOLAR SYSTEM BODY.
*
*          POS1   = POSITION VECTOR TOWARD OBSERVED OBJECT, WITH RESPECT
*                   TO ORIGIN AT OBSERVER (OR THE GEOCENTER) (IN)
*          PE     = POSITION VECTOR OF SOLAR SYSTEM BODY, WITH RESPECT
*                   TO ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (IN)
*          DIFLT  = LIGHT TIME TO POINT ON LINE DEFINED BY POS1 THAT IS
*                   CLOSEST TO SOLAR SYSTEM BODY (POSITIVE IF LIGHT
*                   PASSES BODY BEFORE HITTING OBSERVER, I.E., IF
*                   POS1 IS WITHIN 90 DEGREES OF PE)(OUT)
*
*
      DOUBLE PRECISION POS1,PE,DIFLT,C,DIS,U1,DSQRT
      DIMENSION POS1(3), PE(3), U1(3)
      SAVE
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ('C(AU/DAY)',1.D0,C)
      END IF
 
*     FROM POS1, FORM UNIT VECTOR U1 IN DIRECTION OF STAR OR
*     LIGHT SOURCE
      DIS = DSQRT ( POS1(1)**2 + POS1(2)**2 + POS1(3)**2 )
      DO J=1,3
   20   U1(J) = POS1(J) / DIS
        enddo
*     LIGHT-TIME RETURNED IS THE PROJECTION OF VECTOR PE ONTO THE UNIT
*     VECTOR U1 (FORMED FROM POS1), DIVIDED BY THE SPEED OF LIGHT
      DIFLT = ( PE(1)*U1(1) + PE(2)*U1(2) + PE(3)*U1(3) ) / C
 
      RETURN
 
      END
 
 
 
      SUBROUTINE GRVDEF (TJD,LOC,POS1,POBS,POS2)
*
*     THIS SUBROUTINE COMPUTES THE TOTAL GRAVITATIONAL DEFLECTION OF
*     LIGHT FOR THE OBSERVED OBJECT DUE TO THE MAJOR GRAVITATING BODIES
*     IN THE SOLAR SYSTEM.  THIS SUBROUTINE VALID FOR AN OBSERVED BODY
*     WITHIN THE SOLAR SYSTEM AS WELL AS FOR A STAR.  SEE KLIONER
*     (2003), ASTRONOMICAL JOURNAL 125, 1580-1597, SECTION 6.
*
*          TJD    = TDB JULIAN DATE OF OBSERVATION
*          LOC    = CODE FOR LOCATION OF OBSERVER, DETERMINING
*                   WHETHER THE GRAVITATIONAL DEFLECTION DUE TO THE
*                   EARTH ITSELF IS APPLIED (IN)
*                   SET LOC=0 FOR NO EARTH DEFLECTION (NORMALLY MEANS
*                             OBSERVER IS AT GEOCENTER)
*                   SET LOC=1 TO ADD IN EARTH DEFLECTION (NORMALLY
*                             MEANS OBSERVER IS ON OR ABOVE SURFACE
*                             OF EARTH, INCLUDING EARTH ORBIT)
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
*                   TO ICRS AXES, COMPONENTS IN AU (IN)
*          POBS   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   REFERRED TO ICRS AXES, COMPONENTS IN AU (IN)
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
*                   TO ICRS AXES, CORRECTED FOR GRAVITATIONAL
*                   DEFLECTION, COMPONENTS IN AU (OUT)
*
*
      DOUBLE PRECISION TJD,POS1,POBS,POS2,C,RMASS,RMASSE,PBODY,VBODY,
     .     PBODYO,X,TLT,DLT,TCLOSE,DSQRT
      CHARACTER*3 NAME
      DIMENSION POS1(3), POBS(3), POS2(3), NAME(10), ID(10), RMASS(10),
     .     PBODY(3), VBODY(3), PBODYO(3)
      SAVE
 
*     THE FOLLOWING LIST OF NAMES IDENTIFIES WHICH GRAVITATING BODIES
*     (ASIDE FROM THE EARTH) ARE POTENTIALLY USED -- LIST IS TAKEN FROM
*     KLIONER'S TABLE 1, THE ORDER BASED ON AREA OF SKY AFFECTED (COL 2)
      DATA NAME / 'SUN', 'JUP', 'SAT', 'MOO', 'VEN', 'URA', 'NEP',
     .     3*'   ' /
*     CHANGE VALUE OF NBODY TO INCLUDE OR EXCLUDE GRAVITATING BODIES
*     (NBODY=0 MEANS NO DEFLECTION CALCULATED, NBODY=1 MEANS SUN ONLY,
*     NBODY=2 MEANS SUN + JUPITER, ETC.)
      DATA NBODY / 3 /
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ( 'C(AU/DAY)', 1.D0, C )
*         GET ID NUMBERS AND RECIPROCAL MASSES OF GRAVITATING BODIES
          DO 20 I = 1, NBODY
              ID(I) = IDSS ( NAME(I) )
              CALL ASTCON ( 'MASS_'//NAME(I), 1.D0,   RMASS(I) )
  20      CONTINUE
          IDE = IDSS ( 'EARTH' )
          CALL ASTCON ( 'MASS_EARTH', 1.D0,   RMASSE )
      END IF
 
*     INITIALIZE OUTPUT VECTOR OF OBSERVED OBJECT TO EQUAL INPUT VECTOR
      DO 30 J = 1, 3
          POS2(J) = POS1(J)
  30  CONTINUE
*     OPTION FOR NO DEFLECTION
      IF ( NBODY .LE. 0 ) GO TO 50
 
*     COMPUTE LIGHT-TIME TO OBSERVED OBJECT
      TLT = DSQRT ( POS1(1)**2 + POS1(2)**2 + POS1(3)**2 ) / C
 
*     CYCLE THROUGH GRAVITATING BODIES
      DO 40 I = 1, NBODY
 
          IF ( ID(I) .EQ. -9999 ) GO TO 40
 
*         GET POSITION OF GRAVITATING BODY WRT SS BARYCENTER AT TIME TJD
          CALL SOLSYS ( TJD, ID(I), 0,   PBODY, VBODY, IERR )
 
*         GET POSITION OF GRAVITATING BODY WRT OBSERVER AT TIME TJD
          CALL GEOCEN ( PBODY, POBS,   PBODYO, X )
 
*         COMPUTE LIGHT-TIME FROM POINT ON INCOMING LIGHT RAY THAT
*         IS CLOSEST TO GRAVITATING BODY
          CALL DLIGHT ( POS2, PBODYO,   DLT )
 
*         GET POSITION OF GRAVITATING BODY WRT SS BARYCENTER AT TIME
*         WHEN INCOMING PHOTONS WERE CLOSEST TO IT
          TCLOSE = TJD
          IF ( DLT .GT. 0.D0 ) TCLOSE = TJD - DLT
          IF ( TLT .LT. DLT  ) TCLOSE = TJD - TLT
          CALL SOLSYS ( TCLOSE, ID(I), 0,   PBODY, VBODY, IERR )
 
*         COMPUTE DEFLECTION DUE TO GRAVITATING BODY
          CALL GRVD ( POS2, POBS, PBODY, RMASS(I),   POS2 )
 
  40  CONTINUE
 
*     IF OBSERVER IS NOT AT GEOCENTER, ADD IN DEFLECTION DUE TO EARTH
      IF ( LOC .NE. 0 ) THEN
 
*         GET POSITION OF EARTH WRT SS BARYCENTER AT TIME TJD
          CALL SOLSYS ( TJD, IDE, 0,   PBODY, VBODY, IERR )
 
*         COMPUTE DEFLECTION DUE TO EARTH
          CALL GRVD ( POS2, POBS, PBODY, RMASSE,   POS2 )
 
      END IF
 
  50  RETURN
 
      END
 
 
 
      SUBROUTINE GRVD (POS1,POBS,PBODY,RMASS,POS2)
*
*     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR THE DEFLECTION
*     OF LIGHT IN THE GRAVITATIONAL FIELD OF AN ARBITRARY BODY.  ADAPTED
*     FROM MURRAY (1981) MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
*     SEE ALSO FORMULAE IN SECTION B OF THE ASTRONOMICAL ALMANAC, OR
*     KAPLAN ET AL. (1989) ASTRONOMICAL JOURNAL 97, 1197-1210, SECTION
*     III F.  THIS SUBROUTINE VALID FOR AN OBSERVED BODY WITHIN THE
*     SOLAR SYSTEM AS WELL AS FOR A STAR.
*
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (IN)
*          POBS   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU (IN)
*          PBODY  = POSITION VECTOR OF GRAVITATING BODY, WITH RESPECT TO
*                   ORIGIN AT SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU (IN)
*          RMASS  = RECIPROCAL MASS OF GRAVITATING BODY IN SOLAR MASS
*                   UNITS, THAT IS, SUN MASS / BODY MASS (IN)
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), CORRECTED FOR
*                   GRAVITATIONAL DEFLECTION, COMPONENTS IN AU (OUT)
*
*
      DOUBLE PRECISION POS1,POBS,PBODY,RMASS,POS2,C,MAU,GS,PQ,PE,
     .     PMAG,EMAG,QMAG,PHAT,EHAT,QHAT,PDOTQ,EDOTP,QDOTE,
     .     FAC1,FAC2,P2J,DABS,DSQRT
      DIMENSION POS1(3), POBS(3), PBODY(3), POS2(3), PQ(3), PE(3),
     .     PHAT(3), EHAT(3), QHAT(3)
      SAVE
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN METERS/SECOND
          CALL ASTCON ( 'C', 1.D0, C )
*         GET MAU, THE LENGTH OF THE AU IN METERS
          CALL ASTCON ( 'AU', 1.D0, MAU )
*         GET GS, THE HELIOCENTRIC GRAVITATIONAL CONSTANT
          CALL ASTCON ( 'GS', 1.D0, GS )
      END IF
 
*     CONSTRUCT VECTOR PQ FROM GRAVITATING BODY TO OBSERVED OBJECT AND
*     CONSTRUCT VECTOR PE FROM GRAVITATING BODY TO OBSERVER
      DO 20 J=1,3
          PQ(J) = POBS(J) + POS1(J) - PBODY(J)
          PE(J) = POBS(J) - PBODY(J)
  20  CONTINUE
 
*     COMPUTE VECTOR MAGNITUDES AND UNIT VECTORS
      PMAG = DSQRT (POS1(1)**2 + POS1(2)**2 + POS1(3)**2)
      EMAG = DSQRT (  PE(1)**2 +   PE(2)**2 +   PE(3)**2)
      QMAG = DSQRT (  PQ(1)**2 +   PQ(2)**2 +   PQ(3)**2)
      DO 30 J = 1, 3
          PHAT(J) = POS1(J) / PMAG
          EHAT(J) =   PE(J) / EMAG
          QHAT(J) =   PQ(J) / QMAG
  30  CONTINUE
 
*     COMPUTE DOT PRODUCTS OF VECTORS
      PDOTQ = PHAT(1)*QHAT(1) + PHAT(2)*QHAT(2) + PHAT(3)*QHAT(3)
      EDOTP = EHAT(1)*PHAT(1) + EHAT(2)*PHAT(2) + EHAT(3)*PHAT(3)
      QDOTE = QHAT(1)*EHAT(1) + QHAT(2)*EHAT(2) + QHAT(3)*EHAT(3)
 
*     IF GRAVITATING BODY IS OBSERVED OBJECT, OR IS ON A STRAIGHT LINE
*     TOWARD OR AWAY FROM OBSERVED OBJECT TO WITHIN 1 ARCSEC,
*     DEFLECTION IS SET TO ZERO
      IF ( DABS ( EDOTP ) .GT. 0.99999999999D0 ) THEN
          DO 35 J=1,3
              POS2(J) = POS1(J)
  35      CONTINUE
          GO TO 50
      END IF
 
*     COMPUTE SCALAR FACTORS
      FAC1 = 2.0D0 * GS / (C * C * EMAG * MAU * RMASS)
      FAC2 = 1.0D0 + QDOTE
 
*     CONSTRUCT CORRECTED POSITION VECTOR POS2
      DO 40 J = 1, 3
          P2J = PHAT(J) + FAC1 * (PDOTQ*EHAT(J) - EDOTP*QHAT(J)) / FAC2
          POS2(J) = P2J * PMAG
  40  CONTINUE
 
  50  RETURN
 
      END
 
 
 
      SUBROUTINE ABERAT (POS1,VE,TLIGHT,POS2)
*
*     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR ABERRATION OF LIGHT.
*     ALGORITHM INCLUDES RELATIVISTIC TERMS.  ADAPTED FROM MURRAY (1981)
*     MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
*
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH REESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (IN)
*          VE     = VELOCITY VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU/DAY (IN)
*          TLIGHT = LIGHT TIME FROM BODY TO OBSERVER (OR THE GEOCENTER)
*                   IN DAYS (IN)
*                   IF TLIGHT = 0.D0, THIS SUBROUTINE WILL COMPUTE
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), CORRECTED
*                   FOR ABERRATION, COMPONENTS IN AU (OUT)
*
*
      DOUBLE PRECISION POS1,VE,TLIGHT,POS2,C,TL,P1MAG,VEMAG,
     .     BETA,DOT,COSD,GAMMAI,P,Q,R,DSQRT
      DIMENSION POS1(3), VE(3), POS2(3)
      SAVE
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ('C(AU/DAY)',1.D0,C)
      END IF
 
      TL = TLIGHT
      P1MAG = TL * C
      IF (TL.NE.0.D0) GO TO 20
      P1MAG = DSQRT(POS1(1)**2 + POS1(2)**2 + POS1(3)**2)
      TL = P1MAG / C
   20 VEMAG = DSQRT(VE(1)**2 + VE(2)**2 + VE(3)**2)
      BETA = VEMAG / C
      DOT = POS1(1)*VE(1) + POS1(2)*VE(2) + POS1(3)*VE(3)
      COSD = DOT / (P1MAG * VEMAG)
      GAMMAI = DSQRT(1.D0 - BETA**2)
      P = BETA * COSD
      Q = (1.D0 + P / (1.D0 + GAMMAI)) * TL
      R = 1.D0 + P
 
      DO J=1,3
   30    POS2(J) = (GAMMAI * POS1(J) + Q * VE(J)) / R
         enddo
      RETURN
 
      END
 
 
 
      SUBROUTINE RADVL ( POS, VEL, VELOBS, STAR, DIST,   RV )
*
*     THIS SUBROUTINE PREDICTS THE RADIAL VELOCITY OF THE OBSERVED
*     OBJECT AS IT WOULD BE MEASURED BY SPECTROSCOPIC MEANS.  RADIAL
*     VELOCITY IS HERE DEFINED AS THE RADIAL VELOCITY MEASURE (Z)
*     TIMES THE SPEED OF LIGHT.  FOR A SOLAR SYSTEM BODY, IT APPLIES
*     TO A FICTITIOUS EMITTER AT THE CENTER OF THE OBSERVED OBJECT,
*     ASSUMED MASSLESS (NO GRAVITATIONAL RED SHIFT), AND DOES NOT
*     IN GENERAL APPLY TO REFLECTED LIGHT.  FOR STARS, IT INCLUDES
*     ALL EFFECTS, SUCH AS GRAVITATIONAL RED SHIFT, CONTAINED
*     IN THE CATALOG BARYCENTRIC RADIAL VELOCITY MEASURE, A SCALAR
*     DERIVED FROM SPECTROSCOPY.  NEARBY STARS WITH A KNOWN KINEMATIC
*     VELOCITY VECTOR (OBTAINED INDEPENDENTLY OF SPECTROSCOPY) CAN BE
*     TREATED LIKE SOLAR SYSTEM OBJECTS.  SEE LINDEGREN & DRAVINS
*     (2003), ASTRONOMY & ASTROPHYSICS 401, 1185-1201.
*
*          POS    = GEOMETRIC POSITION VECTOR OF OBJECT WITH RESPECT TO
*                   OBSERVER, CORRECTED FOR LIGHT-TIME, IN AU (IN)
*          VEL    = VELOCITY VECTOR OF OBJECT WITH RESPECT TO SOLAR
*                   SYSTEM BARYCENTER, COMPONENTS IN AU/DAY (IN)
*          VELOBS = VELOCITY VECTOR OF OBSERVER WITH RESPECT TO SOLAR
*                   SYSTEM BARYCENTER, COMPONENTS IN AU/DAY (IN)
*          STAR   = 3-ELEMENT ARRAY OF CATALOG DATA FOR A STAR, TO BE
*                   NON-ZERO IF OBSERVED OBJECT IS A STAR FOR WHICH THE
*                   CATALOG RADIAL VELOCITY IS CONSISTENT WITH
*                   THE IAU DEFINITION OF BARYCENTRIC RADIAL VELOCITY
*                   MEASURE (OTHERWISE ALL ELEMENTS SHOULD BE SET TO
*                   0.D0 EXACTLY) (IN)
*                   STAR(1) = CATALOG RA IN HOURS
*                   STAR(2) = CATALOG DEC IN DEGREES
*                   STAR(3) = Z*C, THE CATALOG BARYCENTRIC RADIAL
*                             VELOCITY MEASURE TIMES THE SPEED OF LIGHT,
*                             IN KILOMETERS/SECOND
*                   ALL THREE DATA ELEMENTS MUST APPLY TO THE SAME
*                   EPOCH (USUALLY J2000.0 = JD 2451545.0 TT)
*          DIST   = 3-ELEMENT ARRAY OF DISTANCES IN AU (IN)
*                   DIST(1) = DISTANCE OF OBSERVER FROM THE GEOCENTER
*                   DIST(2) = DISTANCE OF OBSERVER FROM THE SUN
*                   DIST(3) = DISTANCE OF OBJECT FROM THE SUN
*          RV     = THE OBSERVED RADIAL VELOCITY MEASURE TIMES
*                   THE SPEED OF LIGHT, IN KILOMETERS/SECOND (OUT)
*
*     NOTE 1:  ALL THE INPUT ARGUMENTS ARE BCRS QUANTITIES, EXPRESSED
*     WITH RESPECT TO THE ICRS AXES.  VEL AND VELOBS ARE KINEMATIC
*     VELOCITIES -- DERIVED FROM GEOMETRY OR DYNAMICS, NOT SPECTROSCOPY.
*
*     NOTE 2:  IF ANY ELEMENT OF ARRAY STAR IS NON-ZERO, THE ALGORITHM
*     USED WILL BE CONSISTENT WITH THE IAU DEFINITION OF STELLAR
*     RADIAL VELOCITY, SPECIFICALLY, THE BARYCENTRIC RADIAL VELOCITY
*     MEASURE, WHICH IS DERIVED FROM SPECTROSCOPY.  IN THAT CASE,
*     THE VECTOR VEL CAN BE VERY APPROXIMATE -- OR, FOR DISTANT STARS
*     OR GALAXIES, ZERO -- AS IT WILL BE USED ONLY FOR A SMALL GEOMETRIC
*     CORRECTION THAT IS PROPORTIONAL TO PROPER MOTION.
*
*     NOTE 3:  ANY OF THE DISTANCES IN ARRAY DIST CAN BE SET TO ZERO
*     (0.D0) IF THE CORRESPONDING GENERAL RELATIVISTIC GRAVITATIONAL
*     POTENTIAL TERM IS NOT TO BE EVALUATED.  THESE TERMS
*     GENERALLY ARE IMPORTANT ONLY AT THE METER/SECOND LEVEL.  IF
*     THE FIRST TWO DISTANCES ARE BOTH ZERO, AN AVERAGE VALUE
*     WILL BE USED FOR THE RELATIVISTIC TERM FOR THE OBSERVER,
*     APPROPRIATE FOR AN OBSERVER ON THE SURFACE OF THE EARTH.  THE
*     THIRD DISTANCE, IF GIVEN, IS USED ONLY FOR SOLAR SYSTEM OBJECTS.
*
*
      DOUBLE PRECISION POS,VEL,VELOBS,STAR,DIST,RV,PI,RADCON,
     .     AU,C,GS,GE,C2,TOMS,TOMS2,
     .     POSMAG,UK,V2,VO2,R,PHIGEO,PHISUN,REL,ZC,RA,DC,DU,
     .     ZB1,KVOBS,KV,ZOBS1,DSQRT,DCOS,DSIN
 
      LOGICAL DOSTAR
 
      DIMENSION POS(3), VEL(3), VELOBS(3), STAR(3), DIST(3), UK(3),
     .     DU(3)
 
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
 
      IF ( NTIMES .EQ. 1 ) THEN
*         GET AU, LENGTH OF ASTRONOMICAL UNIT IN METERS
          CALL ASTCON ( 'AU', 1.D0,   AU )
*         GET C, THE SPEED OF LIGHT IN METERS/SECOND
          CALL ASTCON ( 'C', 1.D0,   C )
*         GET GS, HELIOCENTRIC GRAVITATIONAL CONSTANT
          CALL ASTCON ( 'GS', 1.D0,   GS )
*         GET GE, GEOCENTRIC GRAVITATIONAL CONSTANT
          CALL ASTCON ( 'GE', 1.D0,   GE )
*         (GS AND GE ARE IN METERS**3/SECOND**2)
          C2 = C**2
          TOMS = AU / 86400.D0
          TOMS2 = TOMS**2
      END IF
 
      RV = 0.D0
 
*     COMPUTE LENGTH OF POSITION VECTOR = DISTANCE TO OBJECT IN AU
      POSMAG = DSQRT ( POS(1)**2 + POS(2)**2 + POS(3)**2 )
      IF ( POSMAG .LT. 1.D-8 ) GO TO 70
 
*     DETERMINE HOW OBJECT IS TO BE PROCESSED
      DOSTAR = STAR(1) .NE. 0.D0 .OR.
     .         STAR(2) .NE. 0.D0 .OR.
     .         STAR(3) .NE. 0.D0
 
*     COMPUTE UNIT VECTOR TOWARD OBJECT
      DO 40 J = 1, 3
          UK(J) = POS(J) / POSMAG
  40  CONTINUE
 
*     COMPUTE VELOCITY-SQUARED FACTORS
      V2  = ( VEL(1)   **2 + VEL(2)   **2 + VEL(3)   **2 ) * TOMS2
      VO2 = ( VELOBS(1)**2 + VELOBS(2)**2 + VELOBS(3)**2 ) * TOMS2
 
*     COMPUTE GEOPOTENTIAL AT OBSERVER, UNLESS OBSERVER IS GEOCENTRIC
      R = DIST(1) * AU
      PHIGEO = 0.D0
      IF ( R .GT. 1.D6 ) PHIGEO = GE / R
 
*     COMPUTE SOLAR POTENTIAL AT OBSERVER
      R = DIST(2) * AU
      PHISUN = 0.D0
      IF ( R .GT. 1.D8 ) PHISUN = GS / R
 
*     COMPUTE RELATIVISTIC POTENTIAL AND VELOCITY FACTOR FOR OBSERVER
      IF ( DIST(1) .NE. 0.D0 .OR. DIST(2) .NE. 0.D0 ) THEN
*         LINDEGREN & DRAVINS EQ. (41), SECOND FACTOR IN PARENTHESES
          REL = 1.D0 - ( PHIGEO + PHISUN ) / C2 - 0.5D0 * VO2 / C2
      ELSE
*         LINDEGREN & DRAVINS EQ. (42), INVERSE
          REL = 1.D0 - 1.550D-8
      END IF
 
      IF ( DOSTAR ) THEN
 
*         FOR STARS, UPDATE BARYCENTRIC RADIAL VELOCITY MEASURE FOR
*         CHANGE IN VIEW ANGLE
          RA = STAR(1) * 15.D0 * RADCON
          DC = STAR(2) * RADCON
          DU(1) = UK(1) - ( DCOS ( DC ) * DCOS ( RA ) )
          DU(2) = UK(2) - ( DCOS ( DC ) * DSIN ( RA ) )
          DU(3) = UK(3) - ( DSIN ( DC )               )
          ZC = STAR(3) * 1.D3 +
     .       ( VEL(1) * DU(1) + VEL(2) * DU(2) + VEL(3) * DU(3) ) * TOMS
 
*         COMPUTE OBSERVED RADIAL VELOCITY MEASURE OF A STAR (INVERSE OF
*         LINDEGREN & DRAVINS EQ. (41))
          ZB1 = 1.D0 + ZC / C
          KVOBS = ( UK(1) * VELOBS(1)
     .            + UK(2) * VELOBS(2)
     .            + UK(3) * VELOBS(3) ) * TOMS
          ZOBS1 = ZB1 * REL / ( 1.D0 + KVOBS / C )
 
      ELSE
 
*         COMPUTE SOLAR POTENTIAL AT OBJECT, IF WITHIN SOLAR SYSTEM
          R = DIST(3) * AU
          PHISUN = 0.D0
          IF ( R .GT. 1.D8 .AND. R .LT. 1.D16 ) PHISUN = GS / R
 
*         COMPUTE OBSERVED RADIAL VELOCITY MEASURE OF A PLANET OR OTHER
*         OBJECT -- INCLUDING A NEARBY STAR -- WHERE KINEMATIC
*         BARYCENTRIC VELOCITY VECTOR IS KNOWN AND GRAVITATIONAL
*         RED SHIFT IS NEGLIGIBLE (LINDEGREN & DRAVINS EQ. (40),
*         APPLIED AS PER S. KLIONER PRIVATE COMMUNICATION (2006))
          KV = ( UK(1) * VEL(1)
     .         + UK(2) * VEL(2)
     .         + UK(3) * VEL(3) ) * TOMS
          ZB1 = ( 1.D0 + KV / C ) /
     .          ( 1.D0 - PHISUN / C2 - 0.5D0 * V2 / C2 )
          KVOBS = ( UK(1) * VELOBS(1)
     .            + UK(2) * VELOBS(2)
     .            + UK(3) * VELOBS(3) ) * TOMS
          ZOBS1 = ZB1 * REL / ( 1.D0 + KVOBS / C )
 
      END IF
 
*     CONVERT OBSERVED RADIAL VELOCITY MEASURE TO KILOMETERS/SECOND
      RV = ( ZOBS1 - 1.D0 ) * C / 1000.D0
 
  70  RETURN
 
      END
 
 
 
      SUBROUTINE PRECES (TJD1,POS1,TJD2,POS2)
*
*     THIS SUBROUTINE PRECESSES EQUATORIAL RECTANGULAR COORDINATES FROM
*     ONE EPOCH TO ANOTHER.  THE COORDINATES ARE REFERRED TO THE MEAN
*     DYNAMICAL EQUATOR AND EQUINOX OF THE TWO RESPECTIVE EPOCHS.  SEE
*     EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, PP. 103-104,
*     AND CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412,
*     567-586.
*
*          TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
*          POS1 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                 EQUINOX OF FIRST EPOCH (IN)
*          TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
*          POS2 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                 EQUINOX OF SECOND EPOCH (OUT)
*
*     NOTE:  EITHER TJD1 OR TJD2 MUST BE 2451545.0 (J2000.0) TDB.
*
*
      DOUBLE PRECISION TJD1,TJD2,POS1,POS2,PI,SECCON,T0,TLAST,T,
     .     EPS0,PSIA,OMEGAA,CHIA,SA,CA,SB,CB,SC,CC,SD,CD,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /
 
*     INITIALIZE PRECESSION ROTATION MATRIX AS IDENTITY MATRIX
      DATA XX, XY, XZ / 1.D0, 0.D0, 0.D0 /
      DATA YX, YY, YZ / 0.D0, 1.D0, 0.D0 /
      DATA ZX, ZY, ZZ / 0.D0, 0.D0, 1.D0 /
 
   3  FORMAT ( ' PRECES ERROR: PRECESSION FROM JD ', F10.1, ' TO ',
     .     F10.1, ' NOT TO/FROM J2000' )
 
      IF ( TJD1 .NE. T0 .AND. TJD2 .NE. T0 ) THEN
          WRITE ( *, 3 ) TJD1, TJD2
          GO TO 50
      END IF
 
*     T IS TIME IN TDB CENTURIES BETWEEN THE TWO EPOCHS
      T = ( TJD2 - TJD1 ) / 36525.D0
      IF ( TJD2 .EQ. T0 ) T = -T
      IF ( DABS ( T - TLAST ) .LT. 1.D-15 ) GO TO 20
 
*     NUMERICAL COEFFICIENTS OF PSI_A, OMEGA_A, AND CHI_A, ALONG WITH
*     EPSILON_0, THE OBLIQUITY AT J2000.0, ARE 4-ANGLE FORMULATION
*     FROM CAPITAINE ET AL. (2003), EQS. (4), (37), & (39)
      EPS0   = 84381.406D0
      PSIA   = ( ( ( ( -    0.0000000951D0   * T
     .                 +    0.000132851D0  ) * T
     .                 -    0.00114045D0   ) * T
     .                 -    1.0790069D0    ) * T
     .                 + 5038.481507D0     ) * T
      OMEGAA = ( ( ( ( +    0.0000003337D0   * T
     .                 -    0.000000467D0  ) * T
     .                 -    0.00772503D0   ) * T
     .                 +    0.0512623D0    ) * T
     .                 -    0.025754D0     ) * T + EPS0
      CHIA   = ( ( ( ( -    0.0000000560D0   * T
     .                 +    0.000170663D0  ) * T
     .                 -    0.00121197D0   ) * T
     .                 -    2.3814292D0    ) * T
     .                 +   10.556403D0     ) * T
      EPS0 = EPS0 / SECCON
      PSIA = PSIA / SECCON
      OMEGAA = OMEGAA / SECCON
      CHIA = CHIA / SECCON
      SA = DSIN ( EPS0 )
      CA = DCOS ( EPS0 )
      SB = DSIN ( -PSIA )
      CB = DCOS ( -PSIA )
      SC = DSIN ( -OMEGAA )
      CC = DCOS ( -OMEGAA )
      SD = DSIN ( CHIA )
      CD = DCOS ( CHIA )
 
*     COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
*     EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)
      XX =  CD * CB - SB * SD * CC
      YX =  CD * SB * CA + SD * CC * CB * CA - SA * SD * SC
      ZX =  CD * SB * SA + SD * CC * CB * SA + CA * SD * SC
      XY = -SD * CB - SB * CD * CC
      YY = -SD * SB * CA + CD * CC * CB * CA - SA * CD * SC
      ZY = -SD * SB * SA + CD * CC * CB * SA + CA * CD * SC
      XZ =  SB * SC
      YZ = -SC * CB * CA - SA * CC
      ZZ = -SC * CB * SA + CC * CA
 
      TLAST = T
 
   20 IF ( TJD2 .EQ. T0 ) GO TO 30
 
*     PERFORM ROTATION FROM J2000.0 TO EPOCH
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50
 
*     PERFORM ROTATION FROM EPOCH TO J2000.0
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)
 
   50 RETURN
 
      END
 
 
 
      SUBROUTINE NUTATE (TJD,POS1,POS2)
*
*     THIS SUBROUTINE NUTATES EQUATORIAL RECTANGULAR COORDINATES FROM
*     THE MEAN DYNAMICAL EQUATOR AND EQUINOX OF EPOCH TO THE TRUE
*     EQUATOR AND EQUINOX OF EPOCH.  SEE EXPLANATORY SUPPLEMENT TO THE
*     ASTRONOMICAL ALMANAC, PP. 114-115.
*
*          TJD    = TDB JULIAN DATE OF EPOCH (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                   EQUINOX OF EPOCH (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF EPOCH (OUT)
*
*     NOTE:  IF TJD IS NEGATIVE, INVERSE NUTATION (TRUE TO MEAN)
*     IS APPLIED.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,TJD1,PI,SECCON,OBLM,OBLT,EQEQ,
     .     DPSI,DEPS,COBM,SOBM,COBT,SOBT,CPSI,SPSI,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
      TJD1 = DABS(TJD)
 
      CALL ETILT ( TJD1,   OBLM, OBLT, EQEQ, DPSI, DEPS )
      OBLM = OBLM * 3600.D0 / SECCON
      OBLT = OBLT * 3600.D0 / SECCON
      DPSI = DPSI / SECCON
      DEPS = DEPS / SECCON
      COBM = DCOS ( OBLM )
      SOBM = DSIN ( OBLM )
      COBT = DCOS ( OBLT )
      SOBT = DSIN ( OBLT )
      CPSI = DCOS ( DPSI )
      SPSI = DSIN ( DPSI )
 
*     COMPUTE ELEMENTS OF NUTATION ROTATION MATRIX
      XX =  CPSI
      YX = -SPSI * COBM
      ZX = -SPSI * SOBM
      XY =  SPSI * COBT
      YY =  CPSI * COBM * COBT + SOBM * SOBT
      ZY =  CPSI * SOBM * COBT - COBM * SOBT
      XZ =  SPSI * SOBT
      YZ =  CPSI * COBM * SOBT - SOBM * COBT
      ZZ =  CPSI * SOBM * SOBT + COBM * COBT
   10 IF ( TJD .LT. 0.D0 ) GO TO 30
 
*     PERFORM ROTATION FROM MEAN TO TRUE
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50
 
*     PERFORM ROTATION FROM TRUE TO MEAN
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)
 
   50 RETURN
 
      END
 
 
 
      SUBROUTINE SPIN (ANGL,POS1,POS2)
*
*     THIS SUBROUTINE TRANSFORMS A VECTOR FROM ONE COORDINATE SYSTEM
*     TO ANOTHER WITH SAME ORIGIN AND AXES ROTATED ABOUT THE
*     Z AXIS.
*
*          ANGL   = ANGLE OF COORDINATE SYSTEM ROTATION, POSITIVE
*                   COUNTERCLOCKWISE WHEN VIEWED FROM +Z,
*                   IN DEGREES (IN)
*          POS1   = POSITION VECTOR (IN)
*          POS2   = POSITION VECTOR EXPRESSED IN NEW COORDINATE
*                   SYSTEM ROTATED ABOUT Z BY ANGLE ANG (OUT)
*
*
      DOUBLE PRECISION ANGL,POS1,POS2,PI,ALAST,ANG,COSANG,SINANG,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
 
      DATA ALAST / -999.D0 /
 
      IF ( DABS ( ANGL - ALAST ) .GT. 1.D-12 ) THEN
 
          ANG = ANGL / 180.D0 * PI
          COSANG = DCOS ( ANG )
          SINANG = DSIN ( ANG )
 
*         ROTATION MATRIX FOLLOWS
          XX =  COSANG
          YX =  SINANG
          ZX =  0.D0
          XY = -SINANG
          YY =  COSANG
          ZY =  0.D0
          XZ =  0.D0
          YZ =  0.D0
          ZZ =  1.D0
 
          ALAST = ANGL
 
      END IF
 
*     PERFORM ROTATION
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
 
      RETURN
 
      END
 
 
 
      SUBROUTINE WOBBLE (TJD,XP,YP,POS1,POS2)
*
*     THIS SUBROUTINE CORRECTS A VECTOR IN THE ITRS (A ROTATING EARTH-
*     FIXED SYSTEM) FOR POLAR MOTION, AND ALSO CORRECTS THE LONGITUDE
*     ORIGIN (BY A TINY AMOUNT) TO THE TERRESTRIAL INTERMEDIATE ORIGIN
*     (TIO).  THE ITRS VECTOR IS THEREBY TRANSFORMED TO THE TERRESTRIAL
*     INTERMEDIATE SYSTEM, BASED ON THE TRUE (ROTATIONAL) EQUATOR AND
*     THE TERRESTRIAL INTERMEDIATE ORIGIN (TIO).  SINCE THE TRUE EQUATOR
*     IS THE PLANE ORTHOGONAL TO THE DIRECTION OF THE CELESTIAL
*     INTERMEDIATE POLE (CIP), THE COMPONENTS OF THE OUTPUT VECTOR ARE
*     REFERRED TO Z AND X AXES TOWARD THE CIP AND TIO, RESPECTIVELY.
*
*          TJD    = TT OR UT1 JULIAN DATE (IN)
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND TIO (OUT)
*
*     NOTE 1:  IF TJD IS NEGATIVE, THE INVERSE TRANSFORMATION (TERRESTRIAL
*     INTERMEDIATE SYSTEM TO ITRS) IS APPLIED.
*
*     NOTE 2: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION TJD,XP,YP,POS1,POS2,PI,SECCON,T0,T,XPOLE,YPOLE,
     .     SPRIME,TIOLON,SINX,COSX,SINY,COSY,SINL,COSL,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DSIN,DCOS
      DIMENSION POS1(3), POS2(3)
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
*     T0 = TT JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.0D0 /
 
      XPOLE = XP / SECCON
      YPOLE = YP / SECCON
 
      T = ( DABS(TJD) - T0 ) / 36525.D0
 
*     COMPUTE APPROXIMATE LONGITUDE OF TIO, USING EQ. (10) OF
*     LAMBERT & BIZOUARD (2002), ASTRONOMY AND ASTROPHYSICS 394,
*     317-321
      SPRIME = -47.0D-6 * T
      TIOLON = -SPRIME / SECCON
*     NOTE THAT TIOLON, THE LONGITUDE CORRECTION, IS NEGLIGIBLE FOR
*     MOST ASTRONOMICAL PURPOSES
 
*     COMPUTE ELEMENTS OF ROTATION MATRIX
*     EQUIVALENT TO R3(-S')R2(X)R1(Y) AS PER IERS CONVENTIONS (2003)
      SINX = DSIN ( XPOLE )
      COSX = DCOS ( XPOLE )
      SINY = DSIN ( YPOLE )
      COSY = DCOS ( YPOLE )
      SINL = DSIN ( TIOLON )
      COSL = DCOS ( TIOLON )
      XX =  COSX * COSL
      YX =  SINX * SINY * COSL + COSY * SINL
      ZX = -SINX * COSY * COSL + SINY * SINL
      XY = -COSX * SINL
      YY =  SINX * SINY * SINL + COSY * COSL
      ZY =  SINX * COSY * SINL + SINY * COSL
      XZ =  SINX
      YZ = -COSX * SINY
      ZZ =  COSX * COSY
   10 IF ( TJD .LT. 0.D0 ) GO TO 30
 
*     PERFORM ROTATION FROM ITRS TO TERRESTRIAL INTERMEDIATE SYSTEM
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50
 
*     PERFORM ROTATION FROM TERRESTRIAL INTERMEDIATE SYSTEM TO ITRS
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)
 
   50 RETURN
 
      END
 
 
 
      SUBROUTINE FRAME (POS1,K,POS2)
*
*     THIS SUBROUTINE TRANSFORMS A VECTOR FROM THE DYNAMICAL REFERENCE
*     SYSTEM TO THE INTERNATIONAL CELESTIAL REFERENCE SYSTEM (ICRS),
*     OR VICE VERSA.  THE DYNAMICAL REFERENCE SYSTEM IS BASED ON THE
*     DYNAMICAL MEAN EQUATOR AND EQUINOX OF J2000.0.  THE ICRS IS
*     BASED ON THE SPACE-FIXED ICRS AXES DEFINED BY THE RADIO CATALOG
*     POSITIONS OF SEVERAL HUNDRED EXTRAGALACTIC OBJECTS.  THE ROTATION
*     MATRIX USED HERE IS EQUIVALENT TO THAT GIVEN BY HILTON AND
*     HOHENKERK (2004), ASTRONOMY AND ASTROPHYSICS 413, 765-770,
*     EQ. (6) AND (8).
*
*          POS1   = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                   COORDINATES (IN)
*          K      = DIRECTION OF ROTATION (IN)
*                   SET K < 0 FOR DYNAMICAL TO ICRS
*                   SET K > 0 FOR ICRS TO DYNAMICAL
*          POS2   = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                   COORDINATES (OUT)
*
*     NOTE:  FOR GEOCENTRIC COORDINATES, THE SAME TRANSFORMATION IS
*     USED BETWEEN THE DYNAMICAL REFERENCE SYSTEM AND THE GCRS.
*
*
      DOUBLE PRECISION POS1,POS2,PI,SECCON,XI0,ETA0,DA0,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ
      DIMENSION POS1(3), POS2(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
*     XI0, ETA0, AND DA0 ARE ICRS FRAME BIASES IN ARCSECONDS TAKEN
*     FROM IERS CONVENTIONS (2003), CHAPTER 5
      DATA XI0, ETA0, DA0 / -0.0166170D0, -0.0068192D0, -0.01460D0 /
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
 
*     COMPUTE ELEMENTS OF ROTATION MATRIX (TO FIRST ORDER)
      IF ( NTIMES .GT. 1 ) GO TO 20
      XX =  1.D0
      YX = -DA0  / SECCON
      ZX =  XI0  / SECCON
      XY =  DA0  / SECCON
      YY =  1.D0
      ZY =  ETA0 / SECCON
      XZ = -XI0  / SECCON
      YZ = -ETA0 / SECCON
      ZZ =  1.D0
*     INCLUDE SECOND-ORDER CORRECTIONS TO DIAGONAL ELEMENTS
      XX = 1.D0 - 0.5D0 * ( YX**2 + ZX**2 )
      YY = 1.D0 - 0.5D0 * ( YX**2 + ZY**2 )
      ZZ = 1.D0 - 0.5D0 * ( ZY**2 + ZX**2 )
   20 IF ( K .GE. 0 ) GO TO 30
 
*     PERFORM ROTATION FROM DYNAMICAL SYSTEM TO ICRS
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50
 
*     PERFORM ROTATION FROM ICRS TO DYNAMICAL SYSTEM
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)
 
   50 RETURN
 
      END
 
 
 
      SUBROUTINE TERRA (GLON,GLAT,HT,ST,POS,VEL)
*
*     THIS SUBROUTINE COMPUTES THE POSITION AND VELOCITY VECTORS OF
*     A TERRESTRIAL OBSERVER WITH RESPECT TO THE GEOCENTER.
*
*          GLON   = LONGITUDE OF OBSERVER WITH RESPECT TO REFERENCE
*                   MERIDIAN (EAST +) IN DEGREES (IN)
*          GLAT   = GEODETIC LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          ST     = LOCAL APPARENT SIDEREAL TIME AT REFERENCE MERIDIAN
*                   IN HOURS (IN)
*          POS    = POSITION VECTOR OF OBSERVER WITH RESPECT TO
*                   GEOCENTER, EQUATORIAL RECTANGULAR COORDINATES,
*                   REFERRED TO TRUE EQUATOR AND EQUINOX OF DATE,
*                   COMPONENTS IN AU (OUT)
*          VEL    = VELOCITY VECTOR OF OBSERVER WITH RESPECT TO
*                   GEOCENTER, EQUATORIAL RECTANGULAR COORDINATES,
*                   REFERRED TO TRUE EQUATOR AND EQUINOX OF DATE,
*                   COMPONENTS IN AU/DAY (OUT)
*
*     NOTE 1:  IF REFERENCE MERIDIAN IS GREENWICH AND ST=0.D0, POS
*     IS EFFECTIVELY REFERRED TO EQUATOR AND GREENWICH.
*
*     NOTE 2:  THIS SUBROUTINE IGNORES POLAR MOTION, UNLESS THE
*     OBSERVER'S LONGITUDE AND LATITUDE HAVE BEEN CORRECTED FOR IT,
*     AND VARIATION IN THE LENGTH OF DAY (ANGULAR VELOCITY OF EARTH).
*     NEGLECT OF POLAR MOTION MAY YIELD 15 METERS ERROR IN POSITION
*     AND OF ORDER 1 MILLIMETER/SEC ERROR IN VELOCITY.  NEGLECT OF
*     VARIATIONS IN LENGTH OF DAY RESULTS IN EVEN SMALLER VELOCITY
*     ERRORS.
*
*     NOTE 3:  THE TRUE EQUATOR AND EQUINOX OF DATE DO NOT FORM AN
*     INERTIAL SYSTEM.  THEREFORE, WITH RESPECT TO AN INERTIAL SYSTEM,
*     THE SMALL VELOCITY COMPONENT, OF ORDER 0.1 MILLIMETER/SEC,
*     DUE TO THE PRECESSION AND NUTATION OF THE EARTH'S AXIS, IS NOT
*     ACCOUNTED FOR HERE.
*
*
      DOUBLE PRECISION GLON,GLAT,HT,ST,POS,VEL,PI,SECCON,ERAD,F,OMEGA,
     .     AUKM,DF2,PHI,SINPHI,COSPHI,C,S,ACH,ASH,STLOCL,SINST,COSST,
     .     DSQRT,DCOS,DSIN
      DIMENSION POS(3), VEL(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET ERAD, THE EQUATORIAL RADIUS OF EARTH IN KILOMETERS
          CALL ASTCON ('ERAD',1.D-3,ERAD)
*         GET F, THE FLATTENING FACTOR OF EARTH ELLIPSOID
          CALL ASTCON ('F',1.D0,F)
*         GET OMEGA, THE NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF
*         EARTH IN RADIANS/SECOND
          CALL ASTCON ('ANGVEL',1.D0,OMEGA)
*         GET AUKM, THE LENGTH OF THE ASTRONOMICAL UNIT IN KILOMETERS
          CALL ASTCON ('AU',1.D-3,AUKM)
      END IF
 
*     COMPUTE PARAMETERS RELATING TO GEODETIC TO GEOCENTRIC CONVERSION
      DF2 = (1.D0 - F)**2
      PHI = GLAT * 3600.D0 / SECCON
      SINPHI = DSIN(PHI)
      COSPHI = DCOS(PHI)
      C = 1.D0 / DSQRT ( COSPHI**2 + DF2 * SINPHI**2 )
      S = DF2 * C
      ACH = ERAD * C + HT/1000.D0
      ASH = ERAD * S + HT/1000.D0
 
*     COMPUTE LOCAL SIDEREAL TIME FACTORS
      STLOCL = (ST * 54000.D0 + GLON * 3600.D0) / SECCON
      SINST = DSIN(STLOCL)
      COSST = DCOS(STLOCL)
 
*     COMPUTE POSITION VECTOR COMPONENTS IN KM
      POS(1) = ACH * COSPHI * COSST
      POS(2) = ACH * COSPHI * SINST
      POS(3) = ASH * SINPHI
 
*     COMPUTE VELOCITY VECTOR COMPONENTS IN KM/SEC
      VEL(1) = -OMEGA * ACH * COSPHI * SINST
      VEL(2) =  OMEGA * ACH * COSPHI * COSST
      VEL(3) =  0.D0
 
*     CONVERT POSITION AND VELOCITY COMPONENTS TO AU AND AU/DAY
      DO 20 J=1,3
      POS(J) = POS(J) / AUKM
      VEL(J) = VEL(J) / AUKM * 86400.D0
   20 CONTINUE
 
      RETURN
 
      END
 
 
 
      SUBROUTINE TIMES (TDBJD,TTJD,SECDIF)
*
*     THIS SUBROUTINE COMPUTES THE TERRESTRIAL TIME (TT) JULIAN DATE
*     CORRESPONDING TO A BARYCENTRIC DYNAMICAL TIME (TDB) JULIAN DATE.
*     THE EXPRESSION USED IN THIS VERSION IS A TRUNCATED FORM OF A
*     LONGER AND MORE PRECISE SERIES GIVEN BY FAIRHEAD & BRETAGNON
*     (1990) A&A 229, 240.  THE RESULT IS GOOD TO ABOUT 10 MICROSECONDS.
*
*          TDBJD  = TDB JULIAN DATE (IN)
*          TTJD   = TT JULIAN DATE (OUT)
*          SECDIF = DIFFERENCE TDBJD-TTJD, IN SECONDS (OUT)
*
*
      DOUBLE PRECISION TDBJD,TTJD,SECDIF,T,T0,DSIN
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
 
      T = ( TDBJD - T0 ) / 36525.D0
 
*     EXPRESSION GIVEN IN USNO CIRCULAR 179, EQ. 2.6
      SECDIF = 0.001657D0 * DSIN (  628.3076D0 * T + 6.2401D0)
     .       + 0.000022D0 * DSIN (  575.3385D0 * T + 4.2970D0)
     .       + 0.000014D0 * DSIN ( 1256.6152D0 * T + 6.1969D0)
     .       + 0.000005D0 * DSIN (  606.9777D0 * T + 4.0212D0)
     .       + 0.000005D0 * DSIN (   52.9691D0 * T + 0.4444D0)
     .       + 0.000002D0 * DSIN (   21.3299D0 * T + 5.5431D0)
     .       + 0.000010D0 * T * DSIN ( 628.3076D0 * T + 4.2490D0)
 
      TTJD = TDBJD - SECDIF / 86400.D0
 
      RETURN
 
      END
 
 
 
      SUBROUTINE ETILT (TJD,OBLM,OBLT,EQEQ,DPSI,DEPS)
*
*     THIS SUBROUTINE COMPUTES QUANTITIES RELATED TO THE ORIENTATION
*     OF THE EARTH'S ROTATION AXIS AT JULIAN DATE TJD.
*
*          TJD    = TDB JULIAN DATE FOR ORIENTATION PARAMETERS (IN)
*          OBLM   = MEAN OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
*                   DATE TJD (OUT)
*          OBLT   = TRUE OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
*                   DATE TJD (OUT)
*          EQEQ   = EQUATION OF THE EQUINOXES IN TIME SECONDS AT
*                   DATE TJD (OUT)
*          DPSI   = NUTATION IN LONGITUDE IN ARCSECONDS AT
*                   DATE TJD (OUT)
*          DEPS   = NUTATION IN OBLIQUITY IN ARCSECONDS AT
*                   DATE TJD (OUT)
*
*     NOTE:  THE EQUATION OF THE EQUINOXES INCLUDES THE COMPLEMENTARY
*     TERMS.
*
*
      DOUBLE PRECISION TJD,OBLM,OBLT,EQEQ,DPSI,DEPS,PI,SECCON,
     .     T0,TLAST,T,PSI,EPS,PSICOR,EPSCOR,CTERMS,DELPSI,DELEPS,
     .     EL,ELP,F,D,OMEGA,OBM,OBT,EE,
     .     DPOLE1,DPOLE2,DX,DY,DZ,SINE,X,DP1,DP2,DP3,
     .     OBLIQ,DABS,DSIN,DCOS,EECT2000
      INTEGER ACCDIF
      DIMENSION DP1(3), DP2(3), DP3(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   MLAST / 0 /
      DATA DELPSI, DELEPS, CTERMS, PSICOR, EPSCOR / 5 * 0.D0 /
 
*     FUNCTION TO COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
*     CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412, 567-586,
*     EXPRESSION FROM EQ. (39) WITH OBLIQUITY AT J2000.0 TAKEN FROM
*     EQ. (37) OR TABLE 8
      OBLIQ(T) = ( ( ( ( -  0.0000000434D0   * T
     .                   -  0.000000576D0  ) * T
     .                   +  0.00200340D0   ) * T
     .                   -  0.0001831D0    ) * T
     .                   - 46.836769D0     ) * T + 84381.406D0
 
*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )
 
*     CHECK FOR DIFFERENCE IN ACCURACY MODE FROM LAST CALL
      ACCDIF = MOD ( MODE, 2 ) - MOD ( MLAST, 2 )
 
      T = ( TJD - T0 ) / 36525.D0
 
      IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 .OR. ACCDIF .NE. 0 ) THEN
 
*         OBTAIN NUTATION PARAMETERS IN ARCSECONDS
          CALL NOD ( T,   PSI, EPS )
 
*         OBTAIN COMPLEMENTARY TERMS FOR EQUATION OF THE EQUINOXES
*         IN ARCSECONDS
          IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN
*             HIGH-ACCURACY MODE
              CTERMS = EECT2000 ( TJD, 0.D0 ) * SECCON
          ELSE
*             LOW-ACCURACY MODE
              CALL FUNARG ( T,   EL, ELP, F, D, OMEGA )
*             SERIES FROM IERS CONVENTIONS (2003), CHAPTER 5,
*             TABLE 5.2C, WITH SOME ADJUSTMENTS TO COEFFICIENT VALUES
*             COPIED FROM IERS FUNCTION EECT2000, WHICH HAS A MORE
*             COMPLETE SERIES
              CTERMS =
     .          2640.96D-6 * DSIN ( OMEGA )
     .        +   63.52D-6 * DSIN ( 2.D0 * OMEGA )
     .        +   11.75D-6 * DSIN ( 2.D0 * F - 2.D0 * D + 3.D0 * OMEGA )
     .        +   11.21D-6 * DSIN ( 2.D0 * F - 2.D0 * D +        OMEGA )
     .        -    4.55D-6 * DSIN ( 2.D0 * F - 2.D0 * D + 2.D0 * OMEGA )
     .        +    2.02D-6 * DSIN ( 2.D0 * F            + 3.D0 * OMEGA )
     .        +    1.98D-6 * DSIN ( 2.D0 * F            +        OMEGA )
     .        -    1.72D-6 * DSIN ( 3.D0 * OMEGA )
     .        -    0.87D-6 * T * DSIN ( OMEGA )
*             (TERMS SMALLER THAN 2 MICROARCSECONDS OMITTED)
          END IF
          TLAST = TJD
          MLAST = MODE
 
      END IF
 
      DELPSI = PSI + PSICOR
      DELEPS = EPS + EPSCOR
 
*     COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
      OBM = OBLIQ(T)
 
*     COMPUTE TRUE OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
      OBT = OBM + DELEPS
 
*     COMPUTE EQUATION OF THE EQUINOXES IN ARCSECONDS
      EE = DELPSI * DCOS ( OBM / SECCON ) + CTERMS
 
*     CONVERT TO OUTPUT UNITS
      OBLM = OBM / 3600.D0
      OBLT = OBT / 3600.D0
      EQEQ = EE  / 15.D0
      DPSI = DELPSI
      DEPS = DELEPS
 
      RETURN
 
 
      ENTRY CELPOL (TJD,ITYPE,DPOLE1,DPOLE2)
*
*     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF CELESTIAL POLE
*     OFFSETS FOR HIGH-PRECISION APPLICATIONS.  EACH SET OF OFFSETS IS
*     A CORRECTION TO THE MODELED POSITION OF THE POLE FOR A SPECIFIC
*     DATE, DERIVED FROM OBSERVATIONS AND PUBLISHED BY THE IERS.
*     THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE ANY OTHER ROUTINES
*     FOR A GIVEN DATE.  VALUES OF THE POLE OFFSETS SPECIFIED VIA A CALL
*     TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
*
*          TJD    = TDB OR TT JULIAN DATE FOR POLE OFFSETS (IN)
*          ITYPE  = TYPE OF POLE OFFSET (IN)
*                   SET ITYPE=1 FOR CORRECTIONS TO ANGULAR COORDINATES
*                               OF MODELED POLE REFERRED TO MEAN
*                               ECLIPTIC OF DATE, THAT IS,
*                               DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON
*                   SET ITYPE=2 FOR CORRECTIONS TO COMPONENTS OF
*                               MODELED POLE UNIT VECTOR WITH REFERRED
*                               TO GCRS AXES, THAT IS, DX AND DY
*          DPOLE1 = VALUE OF CELESTIAL POLE OFFSET IN FIRST COORDINATE,
*                   (DELTA-DELTA-PSI OR DX) IN MILLIARCSECONDS (IN)
*          DPOLE2 = VALUE OF CELESTIAL POLE OFFSET IN SECOND COORDINATE,
*                   (DELTA-DELTA-EPSILON OR DY) IN MILLIARCSECONDS (IN)
*
*     NOTE 1:  TJD IS USED ONLY FOR ITYPE=2, TO TRANSFORM DX AND DY TO
*     THE EQUIVALENT DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON VALUES.
*
*     NOTE 2:  FOR ITYPE=2, DX AND DY ARE UNIT VECTOR COMPONENT
*     CORRECTIONS, BUT ARE EXPRESSED IN MILLIARCSECONDS SIMPLY BY
*     MULTIPLYING BY 206264806, THE NUMBER OF MILLIARCSECONDS IN ONE
*     RADIAN.
*
*
      IF ( ITYPE .EQ. 1 ) THEN
 
          PSICOR = DPOLE1 * 1.D-3
          EPSCOR = DPOLE2 * 1.D-3
 
      ELSE
 
          DX = DPOLE1
          DY = DPOLE2
 
          T = ( TJD - T0 ) / 36525.D0
*         COMPUTE SINE OF MEAN OBLIQUITY OF DATE
          SINE = DSIN ( OBLIQ(T) / SECCON )
 
*         THE FOLLOWING ALGORITHM, TO TRANSFORM DX AND DY TO
*         DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON, IS FROM G. KAPLAN
*         (2003), USNO/AA TECHNICAL NOTE 2003-03, EQS. (7)-(9).
 
*         TRIVIAL MODEL OF POLE TRAJECTORY IN GCRS ALLOWS COMPUTATION
*         OF DZ
          X = ( 2004.19D0 * T ) / SECCON
          DZ = - ( X + 0.5D0 * X**3 ) * DX
 
*         FORM POLE OFFSET VECTOR (OBSERVED - MODELED) IN GCRS
          DP1(1) = DX * 1.D-3 / SECCON
          DP1(2) = DY * 1.D-3 / SECCON
          DP1(3) = DZ * 1.D-3 / SECCON
 
*         PRECESS POLE OFFSET VECTOR TO MEAN EQUATOR AND EQUINOX OF DATE
          CALL FRAME ( DP1, 1, DP2 )
          CALL PRECES ( T0, DP2, TJD, DP3 )
 
*         COMPUTE DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON IN ARCSECONDS
          PSICOR = ( DP3(1) / SINE ) * SECCON
          EPSCOR = ( DP3(2)        ) * SECCON
 
      END IF
 
      RETURN
 
      END
 
 
 
      SUBROUTINE FUNARG (T,EL,ELP,F,D,OMEGA)
*
*     THIS SUBROUTINE COMPUTES FUNDAMENTAL ARGUMENTS (MEAN ELEMENTS)
*     OF THE SUN AND MOON.  SEE SIMON ET AL. (1994) ASTRONOMY AND
*     ASTROPHYSICS 282, 663-683, ESPECIALLY SECTIONS 3.4-3.5.
*
*          T      = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          EL     = MEAN ANOMALY OF THE MOON IN RADIANS
*                   AT DATE TJD (OUT)
*          ELP    = MEAN ANOMALY OF THE SUN IN RADIANS
*                   AT DATE TJD (OUT)
*          F      = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
*                   OF THE MOON'S ASCENDING NODE IN RADIANS
*                   AT DATE TJD (OUT)
*          D      = MEAN ELONGATION OF THE MOON FROM THE SUN IN
*                   RADIANS AT DATE TJD (OUT)
*          OMEGA  = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE
*                   IN RADIANS AT DATE TJD (OUT)
*
*
      DOUBLE PRECISION T,EL,ELP,F,D,OMEGA,PI,SECCON,REV,DMOD
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
      PARAMETER ( REV    = 360.D0 * 3600.D0      )
 
*     FUNDAMENTAL (DELAUNAY) ARGUMENTS FROM SIMON ET AL. (1994)
 
*     MEAN ANOMALY OF THE MOON
      EL    = DMOD (         485868.249036D0 +
     .               T*( 1717915923.2178D0 +
     .               T*(         31.8792D0 +
     .               T*(          0.051635D0 +
     .               T*(        - 0.00024470D0 )))), REV ) / SECCON
 
*     MEAN ANOMALY OF THE SUN
      ELP   = DMOD (        1287104.79305D0 +
     .               T*(  129596581.0481D0 +
     .               T*(        - 0.5532D0 +
     .               T*(          0.000136D0 +
     .               T*(        - 0.00001149D0 )))), REV ) / SECCON
 
*     MEAN ARGUMENT OF THE LATITUDE OF THE MOON
      F     = DMOD (         335779.526232D0 +
     .               T*( 1739527262.8478D0 +
     .               T*(       - 12.7512D0 +
     .               T*(       -  0.001037D0 +
     .               T*(          0.00000417D0 )))), REV ) / SECCON
 
*     MEAN ELONGATION OF THE MOON FROM THE SUN
      D     = DMOD (        1072260.70369D0 +
     .               T*( 1602961601.2090D0 +
     .               T*(        - 6.3706D0 +
     .               T*(          0.006593D0 +
     .               T*(        - 0.00003169D0 )))), REV ) / SECCON
 
*     MEAN LONGITUDE OF THE ASCENDING NODE OF THE MOON (FROM SIMON
*     SECTION 3.4(b.3), PRECESSION=5028.8200 ARCSEC/CY)
      OMEGA = DMOD (         450160.398036D0 +
     .               T*(  - 6962890.5431D0 +
     .               T*(          7.4722D0 +
     .               T*(          0.007702D0 +
     .               T*(        - 0.00005939D0 )))), REV ) / SECCON
 
      RETURN
 
      END
 
 
 
      SUBROUTINE REFRAC (HEIGHT,ZDOBS,REFR)
*
*     THIS SUBROUTINE COMPUTES ATMOSPHERIC REFRACTION IN ZENITH
*     DISTANCE.  THIS VERSION COMPUTES APPROXIMATE REFRACTION FOR
*     OPTICAL WAVELENGTHS.  IT CAN BE USED FOR PLANNING OBSERVATIONS
*     OR TELESCOPE POINTING, BUT SHOULD NOT BE USED FOR THE REDUCTION
*     OF PRECISE OBSERVATIONS.  BASIC ALGORITHM IS DESCRIBED IN THE
*     EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, P. 144,
*     AND IS AN ADAPTATION OF A FORMULA IN BENNETT (1982), JOURNAL
*     OF NAVIGATION (ROYAL INSTITUTE) 35, 255-259.
*
*          HEIGHT = HEIGHT OF OBSERVER IN METERS (IN)
*          ZDOBS  = OBSERVED ZENITH DISTANCE IN DEGREES (IN)
*          REFR   = ATMOSPHERIC REFRACTION IN DEGREES (OUT)
*
*     NOTE:  HEIGHT IS NOT USED IF ENTRY REFDAT HAS BEEN CALLED
*     TO SPECIFY ATMOSPHERIC PRESSURE.
*
*
      DOUBLE PRECISION HEIGHT,ZDOBS,REFR,PI,DEGRAD,S,
     .     POBS,TOBS,DOBS,WLOBS,OBSP,OBST,OBSD,OBSWL,P,T,D,WL,H,R,
     .     DEXP,DTAN
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           )
 
*     S IS APPROXIMATE SCALE HEIGHT OF ATMOSPHERE IN METERS
      DATA S / 9.1D3 /
      DATA POBS, TOBS, DOBS, WLOBS / 4 * -999.D0 /
 
*     COMPUTE REFRACTION ONLY FOR ZENITH DISTANCES
*     BETWEEN 0.1 AND 91 DEGREES
      IF ( ZDOBS .LT. 0.1D0 .OR. ZDOBS .GT. 91.D0 ) THEN
          REFR = 0.D0
          GO TO 77
      END IF
 
*     IF OBSERVED WEATHER DATA ARE AVAILABLE, USE THEM
*     OTHERWISE, USE CRUDE ESTIMATES OF AVERAGE CONDITIONS
      IF ( POBS .GE. 1.D0 .AND. TOBS .GT. -100.D0 ) THEN
          P  = POBS
          T  = TOBS
          D  = DOBS
          WL = WLOBS
      ELSE
          P  = 1010.D0 * DEXP ( -HEIGHT / S )
          T  = 10.D0
          D  =  0.D0
          WL =  0.5D0
      END IF
*     D AND WL NOT USED IN THIS VERSION
 
      H = 90.D0 - ZDOBS
      R = 0.016667D0 / DTAN ( ( H +  7.31D0 / ( H + 4.4D0 ) ) * DEGRAD )
      REFR = R * ( 0.28D0 * P / ( T + 273.D0 ) )
 
  77  RETURN
 
 
      ENTRY REFDAT (OBSP,OBST,OBSD,OBSWL)
*
*     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF WEATHER OBSERVATIONS
*     AND OTHER DATA TO BE USED IN THE ATMOSPHERIC REFRACTION
*     CALCULATION.  THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE
*     SUBROUTINE REFRAC OR ZDAZ FOR A GIVEN DATE/TIME.  DATA SPECIFIED
*     VIA A CALL TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
*
*          OBSP   = OBSERVED ATMOSPHERIC PRESSURE IN MILLIBARS (IN)
*          OBST   = OBSERVED TEMPERATURE IN DEGREES CELSIUS (IN)
*          OBSD   = OBSERVED DEW POINT IN DEGREES CELSIUS (IN)
*          OBSWL  = OBSERVING WAVELENGTH IN MICRONS (IN)
*
*     NOTE:  OBSD AND OBSWL ARE NOT USED IN THIS VERSION'S REFRACTION
*     ALGORITHM, AND CAN BE SET TO ANY VALUE.
*
*
      POBS  = OBSP
      TOBS  = OBST
      DOBS  = OBSD
      WLOBS = OBSWL
      RETURN
 
      END
 
 
 
      SUBROUTINE LIMANG (POS1,POSO,ALIMB,AFRAC)
*
*     THIS FUNCTION DETERMINES THE ANGLE OF AN OBJECT ABOVE OR BELOW
*     THE EARTH'S LIMB (HORIZON).  THE GEOMETRIC LIMB IS COMPUTED,
*     ASSUMING THE EARTH TO BE AN AIRLESS SPHERE (NO REFRACTION OR
*     OBLATENESS IS INCLUDED).  THE OBSERVER CAN BE ON OR ABOVE THE
*     EARTH.  FOR AN OBSERVER ON THE SURFACE OF THE EARTH, THIS
*     SUBROUTINE RETURNS THE APPROXIMATE UNREFRACTED ALTITUDE.
*
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT GEOCENTER, COMPONENTS IN AU (IN)
*          POSO   = POSITION VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
*                   AT GEOCENTER, COMPONENTS IN AU (IN)
*          ALIMB  = ANGLE OF OBSERVED OBJECT ABOVE (+) OR BELOW (-) LIMB
*                   IN DEGREES (OUT)
*          AFRAC  = NADIR ANGLE OF OBSERVED OBJECT AS A FRACTION OF
*                   APPARENT RADIUS OF LIMB (OUT)
*                   AFRAC<1.D0 MEANS BELOW THE LIMB
*                   AFRAC=1.D0 MEANS ON THE LIMB
*                   AFRAC>1.D0 MEANS ABOVE THE LIMB
*
*
      DOUBLE PRECISION POS1,POSO,ALIMB,AFRAC,PI,HALFPI,DEGCON,
     .     ERAD,AU,RADE,DISOBJ,DISOBS,APRAD,ZDLIM,COSZD,ZDOBJ,
     .     DSQRT,DACOS
      DIMENSION POS1(3), POSO(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( HALFPI = 0.5D0 * PI            )
      PARAMETER ( DEGCON = 180.D0 / PI           )
 
      DATA NTIMES / 0 /
 
      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
           CALL ASTCON ( 'ERAD', 1.D0,   ERAD )
           CALL ASTCON ( 'AU', 1.D0,   AU )
           RADE = ERAD / AU
      END IF
 
      DISOBJ = DSQRT ( POS1(1)**2 + POS1(2)**2 + POS1(3)**2 )
      DISOBS = DSQRT ( POSO(1)**2 + POSO(2)**2 + POSO(3)**2 )
 
*     COMPUTE APPARENT ANGULAR RADIUS OF EARTH'S LIMB
      IF ( DISOBS .GT. RADE ) THEN
          APRAD = DASIN ( RADE / DISOBS )
      ELSE
          APRAD = HALFPI
      END IF
 
*     COMPUTE ZENITH DISTANCE OF EARTH'S LIMB
      ZDLIM = PI - APRAD
 
*     COMPUTE ZENITH DISTANCE OF OBSERVED OBJECT
      COSZD = ( POS1(1)*POSO(1) + POS1(2)*POSO(2) + POS1(3)*POSO(3) )
     .      / ( DISOBJ * DISOBS )
      IF ( COSZD .LE. -1.D0 ) THEN
         ZDOBJ = PI
      ELSE IF ( COSZD .GE. 1.D0 ) THEN
         ZDOBJ = 0.D0
      ELSE
         ZDOBJ = DACOS ( COSZD )
      END IF
 
*     ANGLE OF OBJECT WRT LIMB IS DIFFERENCE IN ZENITH DISTANCES
      ALIMB = ( ZDLIM - ZDOBJ ) * DEGCON
 
*     NADIR ANGLE OF OBJECT AS A FRACTION OF ANGULAR RADIUS OF LIMB
      AFRAC = ( PI - ZDOBJ ) / APRAD
 
      RETURN
 
      END
 
 
 
      SUBROUTINE CIOLOC ( TJD,   RACIO, K )
*
*     THIS SUBROUTINE RETURNS THE LOCATION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO) FOR A GIVEN JULIAN DATE, AS A
*     RIGHT ASCENSION WITH RESPECT TO EITHER THE GCRS (GEOCENTRIC ICRS)
*     ORIGIN OR THE TRUE EQUINOX OF DATE.  THE CIO IS ALWAYS LOCATED ON
*     THE TRUE EQUATOR (=INTERMEDIATE EQUATOR) OF DATE.
*
*          TJD    = TDB JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (OUT)
*          K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
*                   GIVEN (OUT)
*                   K=1 MEANS GCRS
*                   K=2 MEANS TRUE EQUATOR AND EQUINOX OF DATE
*
*     NOTE:  IF AN EXTERNAL FILE OF CIO RIGHT ASCENSIONS IS AVAILABLE,
*     IT WILL BE USED AND K WILL BE SET TO 1.  OTHERWISE AN INTERNAL
*     COMPUTATION WILL BE USED AND K WILL BE SET TO 2.
*
*
      DOUBLE PRECISION TJD,RACIO,A,TLAST,RLAST,JD,RA,P,EQOR,DABS
      LOGICAL USEFIL
      DIMENSION JD(8), RA(8), A(1)
      SAVE
 
*     NUMBER OF VALUES IN ARRAYS FOR LAGRANGIAN INTERPOLATION
      DATA M / 6 /
 
      DATA TLAST, RLAST, KLAST / 0.D0, 0.D0, 0 /
 
   3  FORMAT ( ' CIOLOC ERROR: CANNOT RETURN CIO RA VALUE FOR JD ',
     .     F10.1 )
 
*     CHECK IF EXTERNAL FILE EXISTS
      CALL CIORD ( 0.D0, 1,   A, A, IERR )
      USEFIL = IERR .EQ. 0
 
*     CHECK IF PREVIOUSLY COMPUTED RA VALUE CAN BE USED
      IF ( DABS ( TJD - TLAST ) .LE. 1.D-8 ) THEN
          RACIO = RLAST
          K = KLAST
          GO TO 77
      END IF
 
* --- IF EXTERNAL FILE EXISTS, INTERPOLATE RA VALUES FROM FILE --------
 
      IF ( USEFIL ) THEN
 
          K = 1
 
*         GET ARRAYS OF VALUES TO INTERPOLATE
          CALL CIORD ( TJD, M,   JD, RA, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) TJD
              RACIO = 99.D0
              GO TO 77
          END IF
 
*         PERFORM LAGRANGIAN INTERPOLATION FOR RA AT TJD
          RACIO = 0.D0
          DO 40 J = 1, M
              P = 1.D0
              DO 30 I = 1, M
                  IF ( I .EQ. J ) GO TO 30
                  P = P * ( TJD - JD(I) ) / ( JD(J) - JD(I) )
  30          CONTINUE
              RACIO = RACIO + P * RA(J)
  40      CONTINUE
 
          RACIO = RACIO / 54000.D0
 
* --- OTHERWISE, USE INTERNAL COMPUTATION ----------------------------
 
      ELSE
 
          K = 2
 
*         GET EQUATION OF THE ORIGINS
          CALL EQXRA ( TJD, 1,   EQOR )
 
          RACIO = -EQOR
 
      END IF
 
* ---------------------------------------------------------------------
 
      TLAST = TJD
      RLAST = RACIO
      KLAST = K
 
  77  RETURN
 
      END
 
 
 
      SUBROUTINE CIORD ( TJD, NVALS,   TLIST, RALIST, IERR )
*
*     GIVEN AN INPUT TDB JULIAN DATE AND THE NUMBER OF DATA POINTS
*     DESIRED, THIS SUBROUTINE RETURNS A SET OF JULIAN DATES AND
*     CORRESPONDING VALUES OF THE GCRS RIGHT ASCENSION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO).  THE RANGE OF DATES IS CENTERED (AT LEAST
*     APPROXIMATELY) ON THE REQUESTED DATE.  THE SUBROUTINE OBTAINS THE
*     DATA FROM AN EXTERNAL DATA FILE.
*
*         TJD    = TDB JULIAN DATE (IN)
*         NVALS  = NUMBER OF JULIAN DATES AND RIGHT ASCENSION VALUES
*                  REQUESTED (NOT LESS THAN 2 OR MORE THAN 20) (IN)
*         TLIST  = ARRAY OF TDB JULIAN DATES (OUT)
*         RALIST = ARRAY OF GCRS RIGHT ASCENSIONS OF THE CIO, FOR THE
*                  JULIAN DATES IN TLIST, IN ARCSECONDS (OUT)
*         IERR   = ERROR INDICATOR (OUT)
*                  IERR=0 MEANS EVERYTHING OK
*                  IERR=1 MEANS TJD BEFORE FIRST USABLE DATE IN FILE
*                  IERR=2 MEANS TJD AFTER LAST USABLE DATE IN FILE
*                  IERR=3 MEANS BAD VALUE OF NVALS
*                  IERR=4 MEANS EXTERNAL FILE CANNOT BE FOUND
*
*     NOTE:  TJD=0.D0 WITH NVALS=1 INDICATES A SPECIAL CALL JUST TO
*     DETERMINE IF EXTERNAL FILE EXISTS.
*
*
      DOUBLE PRECISION TJD,TLIST,RALIST,T,T1,R,TBEG,TEND,TINT,DIF
      CHARACTER FILNAM*40, FILEID*(*)
      LOGICAL FILEOK
      DIMENSION TLIST(NVALS), RALIST(NVALS), T(20), R(20)
      SAVE
 
*     LOGICAL UNIT NUMBER AND FILE ID OF TIME SERIES OF CIO RA VALUES
      DATA LU, ITYP / 24, 1 /
      DATA FILNAM / 'CIO_RA.TXT                              ' /
 
      DATA NTIMES, TBEG, TEND, FILEOK / 0, 0.D0, 1.D10, .FALSE. /
 
   1  FORMAT ( A )
   2  FORMAT ( ' CIORD ERROR: CANNOT FIND FILE ', A )
   3  FORMAT ( ' CIORD ERROR: REQUESTED JD ', F10.1, 1X, A,
     .      ' ALLOWED JD ', F10.1 )
   4  FORMAT ( F16.6, F24.14 )
 
*     SPECIAL CALL JUST TO DETERMINE IF FILE EXITS
*     (NO PRINTED ERROR MESSAGE IF NOT)
      IF ( TJD .EQ. 0.D0 .AND. NVALS .EQ. 1 ) THEN
          IERR = 4
          IF ( NTIMES .EQ. 0 ) INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
          IF ( FILEOK ) IERR = 0
          GO TO 77
      END IF
 
*     IF EXTERNAL FILE IS ALREADY KNOWN NOT TO EXIST, IMMEDIATELY
*     EXIT WITH IERR=4
      IF ( NTIMES .GT. 0 .AND. .NOT. FILEOK ) THEN
          WRITE ( *, 2 ) FILNAM
          IERR = 4
          GO TO 77
      END IF
 
*     CHECK FOR VALID VALUE OF NVALS
      IF ( NVALS .LT. 2 .OR. NVALS .GT. 20 ) THEN
          IERR = 3
          GO TO 77
      END IF
 
      MIDDL = NVALS / 2
 
*     CHECK THAT REQUESTED JULIAN DATE IS WITHIN RANGE OF FILE
  10  IF ( TJD .LT. TBEG ) THEN
          WRITE ( *, 3 )  TJD, 'BEFORE FIRST', TBEG
          IERR = 1
          GO TO 77
      ELSE IF ( TJD .GT. TEND ) THEN
          WRITE ( *, 3 ) TJD, 'AFTER LAST', TEND
          IERR = 2
          GO TO 77
      END IF
 
      IF ( ITYP .EQ. 1 ) THEN
 
*         -------------------------------------------------------------
*         READ JULIAN DATES AND CIO RA VALUES FROM FORMATTED
*         SEQUENTIAL INPUT FILE
 
*         EACH RECORD OF THE FILE MUST CONTAIN A TDB JULIAN DATE
*         AND A CORRESPONDING CIO RA VALUE (WRT GCRS) IN ARCSECONDS
 
*         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL
*         -------------------------------------------------------------
 
*         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
          NTIMES = NTIMES + 1
          IF ( NTIMES .EQ. 1 ) THEN
              INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
              IF ( .NOT. FILEOK ) THEN
                  WRITE ( *, 2 ) FILNAM
                  IERR = 4
                  GO TO 77
              END IF
              OPEN ( UNIT=LU, FILE=FILNAM, FORM='FORMATTED',
     .               ACCESS='SEQUENTIAL', STATUS='OLD' )
              READ ( UNIT=LU, FMT=1 )
              DO 19 I = 1, NVALS
                  READ ( UNIT=LU, FMT=4, END=50 ) T(I), R(I)
  19          CONTINUE
              TINT = NINT ( ( T(2) - T(1) ) * 1000.D0 ) / 1000.D0
              TBEG = T(MIDDL)
              IF ( TJD .LT. TBEG ) GO TO 10
          END IF
 
*         -------------------------------------------------------------
 
*         FILE READ SEQUENCE
 
  20      DIF = ( TJD - T(MIDDL) ) / TINT
          NDIF = DIF
 
*         BASIC DECISION ON HOW TO READ FILE
          IF ( DIF .GE. -0.1D0 .AND. DIF .LE. 1.1D0 ) THEN
*             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
              GO TO 70
          ELSE IF ( DIF .LT. 0.D0 ) THEN
*             REQUESTED JULIAN DATE BEFORE DATA PREVIOUSLY READ
              IREC = ( T(NVALS) - TBEG ) / TINT
              NBACK = 3 * NVALS
              IF ( -DIF .LE. 2 * NVALS .AND. IREC .GT. NBACK ) GO TO 34
              GO TO 25
          ELSE IF ( NDIF .GT. ( NVALS + 1 ) ) THEN
*             REQUESTED JULIAN DATE FAR AHEAD OF DATA PREVIOUSLY READ
              NSKIP = NDIF - NVALS - 1
              GO TO 30
          ELSE
*             REQUESTED JULIAN DATE A BIT AHEAD OF DATA PREVIOUSLY READ
              GO TO 40
          END IF
 
*         REPOSITION FILE AT BEGINNING
  25      REWIND ( UNIT=LU )
          READ ( UNIT=LU, FMT=1 )
          GO TO 36
 
*         FAST SKIP FORWARD
  30      DO 32 I = 1, NSKIP
              READ ( UNIT=LU, FMT=1, END=50 )
  32      CONTINUE
          GO TO 36
 
*         BACKSPACE FILE
  34      DO 35 I = 1, NBACK
              BACKSPACE ( UNIT=LU )
  35      CONTINUE
 
*         FILL UP ARRAYS
  36      DO 38 I = 1, NVALS
              READ ( UNIT=LU, FMT=4, END=50 ) T(I), R(I)
  38      CONTINUE
          GO TO 20
 
*         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
  40      DO 44 I = 1, NVALS - 1
              T(I) = T(I+1)
              R(I) = R(I+1)
  44      CONTINUE
          READ ( UNIT=LU, FMT=4, END=50 ) T(NVALS), R(NVALS)
 
          GO TO 20
 
*         -------------------------------------------------------------
 
*         END OF FILE ENCOUNTERED
  50      BACKSPACE ( UNIT=LU )
          BACKSPACE ( UNIT=LU )
          READ ( UNIT=LU, FMT=4 ) TEND
          TEND = TEND - ( NVALS - MIDDL - 1 ) * TINT
          WRITE ( *, 3 ) TJD, 'AFTER LAST', TEND
          T(MIDDL) = TEND + TINT
          IERR = 2
          GO TO 77
 
      ELSE IF ( ITYP .EQ. 2 ) THEN
 
*         -------------------------------------------------------------
*         READ JULIAN DATES AND CIO RA VALUES FROM UNFORMATTED
*         DIRECT ACCESS INPUT FILE
 
*         EACH RECORD OF THE FILE (EXCEPT THE FIRST) MUST CONTAIN A
*         TDB JULIAN DATE AND A CORRESPONDING CIO RA VALUE (WRT GCRS)
*         IN ARCSECONDS
 
*         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL
 
*         THE FIRST RECORD OF THE FILE IS SPECIAL AND MUST CONTAIN THE
*         TOTAL NUMBER OF RECORDS IN THE FILE
*         -------------------------------------------------------------
 
*         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
          NTIMES = NTIMES + 1
          IF ( NTIMES .EQ. 1 ) THEN
              INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
              IF ( .NOT. FILEOK ) THEN
                  WRITE ( *, 2 ) FILNAM
                  IERR = 4
                  GO TO 77
              END IF
              OPEN ( UNIT=LU, FILE=FILNAM, FORM='UNFORMATTED',
     .               ACCESS='DIRECT', RECL=16, STATUS='OLD' )
              READ ( UNIT=LU, REC=1 ) NRECS
              DO 59 I = 1, NVALS
                  READ ( UNIT=LU, REC=I+1 ) T(I), R(I)
  59          CONTINUE
              TINT = NINT ( ( T(2) - T(1) ) * 1000.D0 ) / 1000.D0
              TBEG = T(MIDDL)
              TEND = T(1) + ( NRECS - NVALS + MIDDL ) * TINT
              T1 = T(1)
              LREC = 1
              MAXREC = NRECS - NVALS + 1
              IF ( TJD .LT. TBEG .OR. TJD .GT. TEND ) GO TO 10
          END IF
 
*         -------------------------------------------------------------
 
*         FILE READ SEQUENCE
 
  60      DIF = ( TJD - T(MIDDL) ) / TINT
*         IREC IS THE DATA RECORD NUMBER (PHYSICAL RECORD NUMBER - 1)
*         OF THE FIRST RECORD IN THE SEQUENCE OF NVALS RECORDS WITH
*         THE RELEVANT DATA TO BE RETURNED
          IREC = ( ( TJD - T1 ) / TINT ) - MIDDL + 1.9D0
          IF ( IREC .LT. 1      ) IREC = 1
          IF ( IREC .GT. MAXREC ) IREC = MAXREC
 
*         BASIC DECISION ON HOW TO READ FILE
          IF ( DIF .GE. -0.1D0 .AND. DIF .LE. 1.1D0 ) THEN
*             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
              GO TO 70
          ELSE IF ( IREC .GT. LREC .AND. IREC - LREC .LE. MIDDL ) THEN
*             REQUESTED JULIAN DATE JUST AHEAD OF DATA PREVIOUSLY READ
              GO TO 62
          ELSE
*             REQUESTED JULIAN DATE IN DIFFERENT PART OF FILE
              GO TO 66
          END IF
 
*         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
  62      DO 64 I = 1, NVALS - 1
              T(I) = T(I+1)
              R(I) = R(I+1)
  64      CONTINUE
          READ ( UNIT=LU, REC=LREC+NVALS+1 ) T(NVALS), R(NVALS)
          LREC = LREC + 1
          GO TO 60
 
*         GO DIRECTLY TO PROPER RECORD RANGE AND FILL UP ARRAYS
  66      DO 68 I = 1, NVALS
              READ ( UNIT=LU, REC=IREC+I ) T(I), R(I)
  68      CONTINUE
          LREC = IREC
 
*         -------------------------------------------------------------
 
      END IF
 
*     GOT DATA, SO FILL OUTPUT ARRAYS
  70  DO 75 I = 1, NVALS
          TLIST(I) = T(I)
          RALIST(I) = R(I)
  75  CONTINUE
      IERR = 0
 
  77  RETURN
 
 
      ENTRY CIOFIL ( LUNIT, FILEID, ITYPE )
*
*     THIS ENTRY ALLOWS SPECIFICATION OF THE LOGICAL UNIT NUMBER AND
*     FILE IDENTIFIER OF THE INPUT FILE CONTAINING A TIME SERIES OF CIO
*     RA VALUES.
*
*          LUNIT  = LOGIAL UNIT NUMBER TO BE USED FOR FILE (IN)
*          FILEID = FILE ID (IN)
*          ITYPE  = TYPE OF FILE (IN)
*                   SET ITYPE=1 FOR FORMATTED SEQUENTIAL FILE
*                   SET ITYPE=2 FOR UNFORMATTED BINARY FILE
*                   SET ITYPE=0 OR ANYTHING ELSE TO CLOSE THE CURRENT
*                               FILE (LUNIT AND FILEID IGNORED)
*
*     NOTE:  AFTER A CALL TO CIOFIL WITH ITYPE=0, THE ORIGINAL OR A
*     DIFFERENT FILE OF CIO RA VALUES CAN BE ACCESSED BY SUBSEQUENT
*     CALLS TO CIORD, BUT ONLY AFTER ANOTHER CALL TO CIOFIL WITH
*     ITYPE=1 OR 2.
*
*
      IF ( ITYPE .EQ. 1 .OR. ITYPE .EQ. 2 ) THEN
          LU = LUNIT
          FILNAM = FILEID
          ITYP = ITYPE
      ELSE
          CLOSE ( UNIT=LU )
          NTIMES = 0
          TBEG = 0.D0
          TEND = 1.D10
          FILEOK = .FALSE.
      END IF
 
      RETURN
 
      END
 
 
 
      SUBROUTINE CIOBAS ( TJD, RACIO, K,   X, Y, Z )
*
*     THIS SUBROUTINE RETURNS THE ORTHONORMAL BASIS VECTORS, WITH
*     RESPECT TO THE GCRS (GEOCENTRIC ICRS), OF THE CELESTIAL
*     INTERMEDIATE SYSTEM DEFINED BY THE CELESTIAL INTERMEDIATE POLE
*     (CIP) (IN THE Z DIRECTION) AND THE CELESTIAL INTERMEDIATE ORIGIN
*     (CIO) (IN THE X DIRECTION).  A TDB JULIAN DATE AND THE RIGHT
*     ASCENSION OF THE CIO AT THAT DATE IS REQUIRED AS INPUT.  THE
*     RIGHT ASCENSION OF THE CIO CAN BE WITH RESPECT TO EITHER THE
*     GCRS ORIGIN OR THE TRUE EQUINOX OF DATE -- DIFFERENT ALGORITHMS
*     ARE USED IN THE TWO CASES.
*
*          TJD    = TDB JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (IN)
*          K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
*                   EXPRESSED (IN)
*                   SET K=1 FOR GCRS
*                   SET K=2 FOR TRUE EQUATOR AND EQUINOX OF DATE
*          X      = UNIT VECTOR TOWARD THE CIO, EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO THE GCRS (OUT)
*          Y      = UNIT VECTOR TOWARD THE Y-DIRECTION, EQUATORIAL
*                   RECTANGULAR COORDINATES, REFERRED TO THE GCRS (OUT)
*          Z      = UNIT VECTOR TOWARD NORTH CELESTIAL POLE (CIP),
*                   EQUATORIAL RECTANGULAR COORDINATES, REFERRED TO
*                   THE GCRS (OUT)
*
*
      DOUBLE PRECISION TJD,RACIO,X,Y,Z,XX,YY,ZZ,W0,W1,W2,Z0,
     .     PI,RADCON,T0,TLAST,SINRA,COSRA,XMAG,DABS,DSIN,DCOS,DSQRT
      DIMENSION X(3), Y(3), Z(3), XX(3), YY(3), ZZ(3), Z0(3),
     .     W0(3), W1(3), W2(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA Z0 / 0.D0, 0.D0, 1.D0 /,   TLAST / 0.D0 /,   KLAST / 0 /
 
   3  FORMAT ( ' CIOBAS ERROR: INVALID VALUE FOR K FOR JD ',
     .     F10.1 )
 
*     USE LAST-COMPUTED BASIS VECTORS IF POSSIBLE
      IF ( DABS ( TJD - TLAST ) .LE. 1.D-8 .AND. K .EQ. KLAST )
     .   GO TO 60
 
*     COMPUTE UNIT VECTOR Z TOWARD CELESTIAL POLE (CIP)
      CALL NUTATE ( -TJD, Z0,   W1 )
      CALL PRECES ( TJD, W1, T0,   W2 )
      CALL FRAME ( W2, -1,   ZZ )
 
* --- RA OF CIO EXPRESSED IN GCRS -------------------------------------
 
      IF ( K .EQ. 1 ) THEN
 
*         COMPUTE VECTOR X TOWARD CIO IN GCRS
          SINRA = DSIN ( RACIO * 15.D0 * RADCON )
          COSRA = DCOS ( RACIO * 15.D0 * RADCON )
          XX(1) =  ZZ(3) * COSRA
          XX(2) =  ZZ(3) * SINRA
          XX(3) = -ZZ(1) * COSRA - ZZ(2) * SINRA
 
*         NORMALIZE VECTOR X
          XMAG = DSQRT ( XX(1)**2 + XX(2)**2 + XX(3)**2 )
          XX(1) = XX(1) / XMAG
          XX(2) = XX(2) / XMAG
          XX(3) = XX(3) / XMAG
 
*         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
          YY(1) = ZZ(2) * XX(3) - ZZ(3) * XX(2)
          YY(2) = ZZ(3) * XX(1) - ZZ(1) * XX(3)
          YY(3) = ZZ(1) * XX(2) - ZZ(2) * XX(1)
 
* --- RA OF CIO EXPRESSED IN EQUATOR-AND-EQUINOX OF DATE SYSTEM -------
 
      ELSE IF ( K .EQ. 2 ) THEN
 
*         CONSTRUCT UNIT VECTOR TOWARD CIO
*         IN EQUATOR-AND-EQUINOX-OF-DATE SYSTEM
          W0(1) = DCOS ( RACIO * 15.D0 * RADCON )
          W0(2) = DSIN ( RACIO * 15.D0 * RADCON )
          W0(3) = 0.D0
 
*         ROTATE THE VECTOR INTO THE GCRS TO FORM UNIT VECTOR X
          CALL NUTATE ( -TJD, W0,   W1 )
          CALL PRECES ( TJD, W1, T0,   W2 )
          CALL FRAME ( W2, -1,   XX )
 
*         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
          YY(1) = ZZ(2) * XX(3) - ZZ(3) * XX(2)
          YY(2) = ZZ(3) * XX(1) - ZZ(1) * XX(3)
          YY(3) = ZZ(1) * XX(2) - ZZ(2) * XX(1)
 
* ---------------------------------------------------------------------
 
      ELSE
 
          WRITE ( *, 3 ) TJD
          GO TO 77
 
      END IF
 
* ---------------------------------------------------------------------
 
      TLAST = TJD
      KLAST = K
 
  60  DO 66 J = 1, 3
          X(J) = XX(J)
          Y(J) = YY(J)
          Z(J) = ZZ(J)
  66  CONTINUE
 
  77  RETURN
 
      END
 
 
 
      SUBROUTINE EROT (DATE1,DATE2,THETA)
*
*     THIS SUBROUTINE RETURNS THE VALUE OF THE EARTH ROTATION ANGLE
*     (THETA) FOR A GIVEN UT1 JULIAN DATE.  THE EXPRESSION USED IS
*     TAKEN FROM THE NOTE TO IAU RESOLUTION B1.8 OF 2000.
*
*         DATE1  = HIGH-ORDER PART OF UT1 JULIAN DATE (IN)
*         DATE2  = LOW-ORDER PART OF UT1 JULIAN DATE (IN)
*         THETA  = EARTH ROTATION ANGLE IN DEGREES (OUT)
*
*
      DOUBLE PRECISION DATE1, DATE2, THETA, T0, THET1, THET2, THET3,
     .     DMOD
 
      DATA T0 / 2451545.0D0 /
 
*     THE ALGORITHM USED BELOW IS EQUIVALENT TO THE CANNONICAL
*     THETA = 0.7790572732640D0 + 1.00273781191135448D0 * T,
*     WHERE T IS THE TIME IN UT1 DAYS FROM 2000.0 (T=DATE1+DATE2-T0),
*     BUT IT AVOIDS MANY TWO-PI 'WRAPS' THAT DECREASE PRECISION
*     (ADOPTED FROM SOFA ROUTINE IAU_ERA00 BY PAT WALLACE; SEE ALSO
*     EXPRESSION AT TOP OF PAGE 35 OF IERS CONVENTIONS (1996))
 
      THET1 = 0.7790572732640D0 + 0.00273781191135448D0 * ( DATE1 - T0 )
      THET2 =                     0.00273781191135448D0 *   DATE2
      THET3 = DMOD ( DATE1, 1.D0 ) + DMOD ( DATE2, 1.D0 )
      THETA = DMOD ( THET1 + THET2 + THET3, 1.D0 ) * 360.D0
      IF ( THETA .LT. 0.D0 ) THETA = THETA + 360.D0
 
      RETURN
 
      END
 
 
 
      SUBROUTINE EQXRA ( TJD, K,    RAEQ )
*
*     THIS SUBROUTINE COMPUTES THE INTERMEDIATE RIGHT ASCENSION
*     OF THE EQUINOX AT JULIAN DATE TJD, USING AN ANALYTICAL EXPRESSION
*     FOR THE ACCUMULATED PRECESSION IN RIGHT ASCENSION.  FOR THE
*     TRUE EQUINOX THE RESULT IS THE EQUATION OF THE ORIGINS.
*
*          TJD    = TDB JULIAN DATE (IN)
*          K      = EQUINOX SELECTION CODE (IN)
*                   SET K=0 FOR MEAN EQUINOX
*                   SET K=1 FOR TRUE EQUINOX (EQUATION OF THE ORIGINS)
*          RADIF  = INTERMEDIATE RIGHT ASCENSION OF THE EQUINOX,
*                   IN HOURS (+ OR -) (OUT)
*
*
      DOUBLE PRECISION TJD,RAEQ,T0,TLAST,EE,EQEQ,T,A,PRECRA,DABS
      SAVE
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   EE / 0.D0 /
 
      T = ( TJD - T0 ) / 36525.D0	
 
*     FOR THE TRUE EQUINOX, OBTAIN THE EQUATION OF THE EQUINOXES IN
*     TIME SECONDS, WHICH INCLUDES THE 'COMPLIMENTARY TERMS'
      IF ( K .EQ. 1 ) THEN
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( TJD,   A, A, EE, A, A )
              TLAST = TJD
          END IF
          EQEQ = EE
      ELSE
          EQEQ = 0.D0
      END IF
 
*     PRECESSION IN RA IN ARCSECONDS TAKEN FROM CAPITAINE ET AL. (2003),
*     ASTRONOMY AND ASTROPHYSICS 412, 567-586, EQ. (42)
      PRECRA = 0.014506D0 +
     .         ( ( ( ( -    0.0000000368D0   * T
     .                 -    0.000029956D0  ) * T
     .                 -    0.00000044D0   ) * T
     .                 +    1.3915817D0    ) * T
     .                 + 4612.156534D0     ) * T
 
      RAEQ = - ( PRECRA / 15.D0 + EQEQ ) / 3600.D0
 
      RETURN
 
      END
 
 
 
      SUBROUTINE SETDT ( DELT )
*
*     THIS SUBROUTINE ALLOWS FOR THE SPECIFICATION OF THE DELTA-T VALUE
*     (DELTA-T = TT - UT1) TO BE USED IN THE CALCULATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS
*     THESE CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE
*     TIME ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME
*     ARGUMENT FOR THE PRECESSION AND NUTATION COMPONENTS.  THIS
*     SUBROUTINE, IF USED, SHOULD BE CALLED BEFORE ANY SUBROUTINE
*     RELATED TO EARTH ROTATION (E.G., SIDTIM OR TERCEL) FOR A GIVEN
*     DATE.  THE VALUE OF DELTA-T SPECIFIED HERE WILL BE USED UNTIL
*     EXPLICITLY CHANGED.
*
*          DELT   = VALUE OF DELTA-T (TT-UT1) IN SECONDS (IN)
*
*     NOTE 1:  THE COMPUTED VALUE OF SIDEREAL TIME, AND THE EQUIVALENT
*     EARTH ORIENTATION ANGLES, ARE RELATIVELY INSENSITIVE TO THE VALUE
*     OF DELTA-T: UP TO ONLY ~3 MICROARCSECONDS PER SECOND OF DELTA-T.
*     THEREFORE, FOR MANY APPLICATIONS, THIS SUBROUTINE EITHER NEED NOT
*     BE CALLED AT ALL, OR CAN BE CALLED JUST ONCE FOR A WIDE RANGE OF
*     DATES (E.G., A YEAR).  IF THIS CALL IS NOT USED, A DEFAULT
*     VALUE OF DELTA-T OF 64 SECONDS IS USED, WHICH IS APPROPRIATE TO
*     2000.0.
*
*     NOTE 2:  THE INPUT TIME ARGUMENTS TO SIDTIM AND TERCEL (TJDH AND
*     TJDL) ARE EXPRESSED IN UT1 REGARDLESS OF WHETHER THIS CALL IS
*     USED.
*
*
      DOUBLE PRECISION DELTAT, DT, DELT
      SAVE DT
 
*     DEFAULT VALUE OF DELTA-T IN DAYS, EQUIVALENT TO 64 SECONDS,
*     THE APPROXIMATE VALUE AT 2000.0
      DATA DT / 0.00074074D0 /
 
      DT = DELT / 86400.D0
 
      RETURN
 
 
      ENTRY GETDT ( DELTAT )
 
*     THIS ENTRY RETURNS THE CURRENT VALUE OF DELTA-T
*     (DELTA-T = TT - UT1), PREVIOUSLY SET BY THE USER.  THE VALUE
*     RETURNED IS TO BE USED IN THE CALCULATION OF SIDEREAL TIME AND
*     THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS THESE
*     CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE TIME
*     ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME ARGUMENT
*     FOR THE PRECESSION AND NUTATION COMPONENTS.
*
*          DELTAT = VALUE OF DELTA-T (TT-UT1) IN DAYS (OUT)
 
 
      DELTAT = DT
 
      RETURN
 
      END
 
 
 
      SUBROUTINE SETMOD ( MODE )
*
*     THIS SUBROUTINE ALLOWS THE USER TO SPECIFY THE 'MODE' VALUE,
*     WHICH DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (IN)
*                   SET MODE=0 FOR CIO-BASED METHOD, FULL ACCURACY
*                   SET MODE=1 FOR CIO-BASED METHOD, REDUCED ACCURACY
*                   SET MODE=2 FOR EQUINOX-BASED METHOD, FULL ACCURACY
*                   SET MODE=3 FOR EQUINOX-BASED METHOD, REDUCED
*                                  ACCURACY
*
*     NOTE: OTHER ENTRY POINTS ARE PROVIDED TO ALLOW THE METHOD AND
*     ACCURACY TO BE SPECIFIED IN A MORE OBVIOUS WAY:
*     MODE=0 CAN BE SET BY CALL CIOTIO AND CALL HIACC
*     MODE=1 CAN BE SET BY CALL CIOTIO AND CALL LOACC
*     MODE=2 CAN BE SET BY CALL EQINOX AND CALL HIACC
*     MODE=3 CAN BE SET BY CALL EQINOX AND CALL LOACC
*
*
      SAVE IMODE, LMODE
 
      DATA IMODE, LMODE / 2, 2 /
 
      LMODE = IMODE
      IMODE = MODE
 
      RETURN
 
 
      ENTRY CIOTIO
      LMODE = IMODE
      IF ( IMODE .GE. 2 ) IMODE = IMODE - 2
      RETURN
 
 
      ENTRY EQINOX
      LMODE = IMODE
      IF ( IMODE .LE. 1 ) IMODE = IMODE + 2
      RETURN
 
 
      ENTRY LOACC
      LMODE = IMODE
      IF ( MOD ( IMODE, 2 ) .EQ. 0 ) IMODE = IMODE + 1
      RETURN
 
 
      ENTRY HIACC
      LMODE = IMODE
      IF ( MOD ( IMODE, 2 ) .EQ. 1 ) IMODE = IMODE - 1
      RETURN
 
 
      ENTRY RESUME
      IMODE = LMODE
      RETURN
 
 
      ENTRY GETMOD ( MODE )
*
*     THIS SUBROUTINE RETURNS THE 'MODE' VALUE, WHICH
*     DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (OUT)
*                   MODE=0 MEANS CIO-BASED METHOD, FULL ACCURACY
*                   MODE=1 MEANS CIO-BASED METHOD, REDUCED ACCURACY
*                   MODE=2 MEANS EQUINOX-BASED METHOD, FULL ACCURACY
*                   MODE=3 MEANS EQUINOX-BASED METHOD, REDUCED ACCURACY
*
*
      MODE = IMODE
 
      RETURN
 
      END
 
 
 
      SUBROUTINE GETVEC ( UNITV )
*
*     THIS SUBROUTINE ALLOWS THE USER TO RETRIEVE THE LAST COMPUTED
*     POSITION ON THE SKY AS A UNIT VECTOR.
*
*          UNITV  = UNIT VECTOR TOWARD LAST COMPUTED POSITION ON THE
*                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
*                   POSITION (OUT)
*
*
      DOUBLE PRECISION UNITV, P, POS, R, DSQRT
      DIMENSION UNITV(3), P(3), POS(3)
      SAVE P
 
      R = DSQRT ( P(1)**2 + P(2)**2 + P(3)**2 )
 
      DO 20 J = 1, 3
          UNITV(J) = P(J) / R
  20  CONTINUE
 
      RETURN
 
 
      ENTRY SETVEC ( POS )
*
*     THIS ENTRY STORES THE LAST COMPUTED POSITION ON THE SKY.
*
*          POS    = VECTOR TOWARD LAST COMPUTED POSITION ON THE
*                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
*                   POSITION (IN)
*
*
      DO 30 J = 1, 3
          P(J) = POS(J)
  30  CONTINUE
 
      RETURN
 
      END
 
 
 
      SUBROUTINE JULDAT (I,M,K,H,TJD)
*
*     THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND
*     TIME.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT TIME VALUE
*     CAN BE IN ANY UT-LIKE TIME SCALE (UTC, UT1, TT, ETC.) - OUTPUT
*     JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
*     VAN FLANDERN.
*
*          I      = YEAR (IN)
*          M      = MONTH NUMBER (IN)
*          K      = DAY OF MONTH (IN)
*          H      = UT HOURS (IN)
*          TJD    = JULIAN DATE (OUT)
*
*
      DOUBLE PRECISION H,TJD
 
*     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
      JD = K-32075+1461*(I+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12
     .     -3*((I+4900+(M-14)/12)/100)/4
      TJD = JD - 0.5D0 + H/24.D0
 
      RETURN
      END
 
 
 
      SUBROUTINE CALDAT (TJD,I,M,K,H)
*
*     THIS SUBROUTINE COMPUTES CALENDAR DATE AND TIME, GIVEN JULIAN
*     DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE TIME SCALE
*     (UTC, UT1, TT, ETC.) - OUTPUT TIME VALUE WILL HAVE SAME BASIS.
*     OUTPUT CALENDAR DATE WILL BE GREGORIAN.  ALGORITHM BY FLIEGEL AND
*     VAN FLANDERN.
*
*          TJD    = JULIAN DATE (IN)
*          I      = YEAR (OUT)
*          M      = MONTH NUMBER (OUT)
*          K      = DAY OF MONTH (OUT)
*          H      = UT HOURS (OUT)
*
*
      DOUBLE PRECISION TJD,H,DJD,DMOD
 
      DJD = TJD + 0.5D0
      JD = DJD
      H = DMOD (DJD,1.D0) * 24.D0
*     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
      L = JD + 68569
      N = 4*L/146097
      L = L - (146097*N+3)/4
*     I=YEAR, M=MONTH, K=DAY
      I = 4000*(L+1)/1461001
      L = L - 1461*I/4 + 31
      M = 80*L/2447
      K = L - 2447*M/80
      L = M / 11
      M = M + 2 - 12*L
      I = 100*(N-49) + I + L
 
      RETURN
      END
 
 
 
      SUBROUTINE ASTCON (NAME,FACTOR,CONST)
*
*     THIS SUBROUTINE SUPPLIES THE VALUES OF ASTRONOMICAL CONSTANTS.
*
*         NAME   = NAME OF CONSTANT WHOSE VALUE IS DESIRED (IN)
*                  'C'         SPEED OF LIGHT IN METERS/SECOND
*                  'C(AU/DAY)' SPEED OF LIGHT IN AU/DAY
*                  'AU'        LENGTH OF ASTRONOMICAL UNIT IN METERS
*                  'AU(SEC)'   LENGTH OF ASTRONOMICAL UNIT IN SECONDS
*                  'GS'        HELIOCENTRIC GRAVITATIONAL CONSTANT
*                                 IN METERS**3/SECOND**2
*                  'GE'        GEOCENTRIC GRAVITATIONAL CONSTANT
*                                 IN METERS**3/SECOND**2
*                  'ERAD'      EQUATORIAL RADIUS OF EARTH IN METERS
*                  'F'         FLATTENING FACTOR OF EARTH
*                  'ANGVEL'    NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY
*                                 OF EARTH IN RADIANS/SECOND
*                  'MASS_SUN'  RECIPROCAL MASS OF THE SUN
*                  'MASS_EAR'  RECIPROCAL MASS OF THE EARTH
*                  'MASS_MOO'  RECIPROCAL MASS OF THE MOON
*                  'MASS_MER'  RECIPROCAL MASS OF MERCURY
*                      :             :      :        :
*                  'MASS_PLU'  RECIPROCAL MASS OF PLUTO
*         FACTOR = FACTOR BY WHICH CONSTANT VALUE IS TO BE MULTIPLIED
*                  (IN)
*         CONST  = CONSTANT VALUE AFTER MULTIPLICATION BY FACTOR (OUT)
*
*
      DOUBLE PRECISION FACTOR,CONST,C,AUSEC
      CHARACTER NAME*(*)
 
*     NOTE:  THESE CONSTANT VALUES ARE BASED ON THE TDB SECOND WHERE
*     APPLICABLE.
 
*     SPEED OF LIGHT IN METERS/SECOND IS A DEFINING PHYSICAL CONSTANT
      DATA C / 299792458.D0 /
 
*     LIGHT-TIME FOR ONE ASTRONOMICAL UNIT IN SECONDS, FROM DE-405
      DATA AUSEC / 499.0047838061D0 /
 
*     SPEED OF LIGHT IN METERS/SECOND
      IF ( NAME .EQ. 'C' ) THEN
          CONST = C
 
*     SPEED OF LIGHT IN AU/DAY
      ELSE IF ( NAME .EQ. 'C(AU/DAY)' ) THEN
          CONST = 86400.D0 / AUSEC
 
*     LENGTH OF ASTRONOMICAL UNIT IN METERS
      ELSE IF ( NAME .EQ. 'AU' ) THEN
          CONST = AUSEC * C
 
*     LENGTH OF ASTRONOMICAL UNIT IN SECONDS
      ELSE IF ( NAME .EQ. 'AU(SEC)' ) THEN
          CONST = AUSEC
 
*     HELIOCENTRIC GRAVITATIONAL CONSTANT IN METERS**3/SECOND**2, FROM
*     DE-405
      ELSE IF ( NAME .EQ. 'GS' ) THEN
          CONST = 1.32712440017987D20
 
*     GEOCENTRIC GRAVITATIONAL CONSTANT IN METERS**3/SECOND**2, FROM
*     DE-405
      ELSE IF ( NAME .EQ. 'GE' ) THEN
          CONST = 3.98600433D14
 
*     EQUATORIAL RADIUS OF EARTH IN METERS, FROM IERS CONVENTIONS (2003)
      ELSE IF ( NAME .EQ. 'ERAD' ) THEN
          CONST = 6378136.6D0
 
*     FLATTENING FACTOR OF EARTH, FROM IERS CONVENTIONS (2003)
      ELSE IF ( NAME .EQ. 'F' ) THEN
          CONST = 1.D0 / 298.25642D0
 
*     NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF EARTH
*     IN RADIANS/SECOND, FROM IERS CONVENTIONS (2003)
      ELSE IF ( NAME .EQ. 'ANGVEL' ) THEN
          CONST = 7.2921150D-5
 
*     RECIPROCAL MASSES OF SOLAR SYSTEM BODIES, FROM DE-405
*     (SUN MASS / BODY MASS)
      ELSE IF ( NAME(1:4) .EQ. 'MASS' ) THEN
 
          CONST = 1.D0
          IF ( NAME(6:8) .EQ. 'SUN' ) CONST =         1.D0
          IF ( NAME(6:8) .EQ. 'MOO' ) CONST =  27068700.387534D0
          IF ( NAME(6:8) .EQ. 'MER' ) CONST =   6023600.D0
          IF ( NAME(6:8) .EQ. 'VEN' ) CONST =    408523.71D0
          IF ( NAME(6:8) .EQ. 'EAR' ) CONST =    332946.050895D0
          IF ( NAME(6:8) .EQ. 'MAR' ) CONST =   3098708.D0
          IF ( NAME(6:8) .EQ. 'JUP' ) CONST =      1047.3486D0
          IF ( NAME(6:8) .EQ. 'SAT' ) CONST =      3497.898D0
          IF ( NAME(6:8) .EQ. 'URA' ) CONST =     22902.98D0
          IF ( NAME(6:8) .EQ. 'NEP' ) CONST =     19412.24D0
          IF ( NAME(6:8) .EQ. 'PLU' ) CONST = 135200000.D0
          IF ( NAME(6:8) .EQ. 'EMB' ) CONST =    328900.561400D0
 
      END IF
 
      CONST = CONST * FACTOR
 
      RETURN
 
      END
 
 
 
      SUBROUTINE NOD (T,DPSI,DEPS)
*
*     THIS SUBROUTINE RETURNS THE VALUES FOR NUTATION IN LONGITUDE AND
*     NUTATION IN OBLIQUITY FOR A GIVEN TDB JULIAN DATE.
*
*          T     = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          DPSI  = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
*          DEPS  = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
*
*
      DOUBLE PRECISION T,DPSI,DEPS,PI,SECCON,T0,T1,DP,DE
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
 
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
 
*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )
 
      T1 = T * 36525.D0
 
*     =================================================================
*     EVALUATE NUTATION SERIES
*     RESULTING NUTATION IN LONGITUDE AND OBLIQUITY IN ARCSECONDS
 
*     CALL SUBROUTINE TO EVALUATE NUTATION SERIES
      IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN
*         HIGH ACCURACY MODE -- IERS SUBROUTINE
          CALL NU2000A ( T0, T1, DP, DE )
      ELSE
*         LOW ACCURACY MODE -- MODIFICATION OF IERS SUBROUTINE
          CALL NU2000K ( T0, T1, DP, DE )
      END IF
      DPSI = DP * SECCON
      DEPS = DE * SECCON
 
*     =================================================================
 
      RETURN
 
      END
 
 
 
      SUBROUTINE NU2000A ( DATE1, DATE2, DPSI, DEPS )
 
*+
*  - - - - - - - -
*   N U 2 0 0 0 A
*  - - - - - - - -
*
*  Nutation, IAU 2000A model (MHB_2000 without FCN).
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
*
*  Returned:
*     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
*
*  This revision:  2002 November 25
*
*-----------------------------------------------------------------------
 
      IMPLICIT NONE
 
      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS
 
*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )
 
*  Milliarcseconds to radians
      DOUBLE PRECISION DMAS2R
      PARAMETER ( DMAS2R = DAS2R / 1D3 )
 
*  Arc seconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )
 
*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )
 
*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )
 
*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )
 
*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )
 
*  Miscellaneous
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 DPSILS, DEPSLS,
     :                 AL, ALSU, AF, AD, AOM, ALME, ALVE, ALEA, ALMA,
     :                 ALJU, ALSA, ALUR, ALNE, APA, DPSIPL, DEPSPL
      INTEGER I, J
 
*  -------------------------
*  Luni-Solar nutation model
*  -------------------------
 
*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 678 )
 
*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)
 
*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)
 
*  ---------------
*  Planetary terms
*  ---------------
 
*  Number of terms in the planetary nutation model
      INTEGER NPL
      PARAMETER ( NPL = 687 )
 
*  Coefficients for fundamental arguments
      INTEGER NAPL(14,NPL)
 
*  Longitude and obliquity coefficients
      INTEGER ICPL(4,NPL)
 
*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------
 
*
*  Luni-Solar argument multipliers:
*               L     L'    F     D     Om
 
      DATA ( ( NALS(I,J), I=1,5 ), J=  1, 10 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 11, 20 ) /
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 21, 30 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 31, 40 ) /
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 41, 50 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 51, 60 ) /
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 61, 70 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 71, 80 ) /
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2,
     :         -2,    0,    2,    4,    2,
     :         -1,    0,    4,    0,    2,
     :          2,    0,    2,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 81, 90 ) /
     :          2,    0,    2,    2,    2,
     :          1,    0,    0,    2,    1,
     :          3,    0,    0,    0,    0,
     :          3,    0,    2,   -2,    2,
     :          0,    0,    4,   -2,    2,
     :          0,    1,    2,    0,    1,
     :          0,    0,   -2,    2,    1,
     :          0,    0,    2,   -2,    3,
     :         -1,    0,    0,    4,    0,
     :          2,    0,   -2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 91,100 ) /
     :         -2,    0,    0,    4,    0,
     :         -1,   -1,    0,    2,    1,
     :         -1,    0,    0,    1,    1,
     :          0,    1,    0,    0,    2,
     :          0,    0,   -2,    0,    1,
     :          0,   -1,    2,    0,    1,
     :          0,    0,    2,   -1,    2,
     :          0,    0,    2,    4,    2,
     :         -2,   -1,    0,    2,    0,
     :          1,    1,    0,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=101,110 ) /
     :         -1,    1,    0,    2,    0,
     :         -1,    1,    0,    1,    2,
     :          1,   -1,    0,    0,    1,
     :          1,   -1,    2,    2,    2,
     :         -1,    1,    2,    2,    2,
     :          3,    0,    2,    0,    1,
     :          0,    1,   -2,    2,    0,
     :         -1,    0,    0,   -2,    1,
     :          0,    1,    2,    2,    2,
     :         -1,   -1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=111,120 ) /
     :          0,   -1,    0,    0,    2,
     :          1,    0,    2,   -4,    1,
     :         -1,    0,   -2,    2,    0,
     :          0,   -1,    2,    2,    1,
     :          2,   -1,    2,    0,    2,
     :          0,    0,    0,    2,    2,
     :          1,   -1,    2,    0,    1,
     :         -1,    1,    2,    0,    2,
     :          0,    1,    0,    2,    0,
     :          0,   -1,   -2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=121,130 ) /
     :          0,    3,    2,   -2,    2,
     :          0,    0,    0,    1,    1,
     :         -1,    0,    2,    2,    0,
     :          2,    1,    2,    0,    2,
     :          1,    1,    0,    0,    1,
     :          1,    1,    2,    0,    1,
     :          2,    0,    0,    2,    0,
     :          1,    0,   -2,    2,    0,
     :         -1,    0,    0,    2,    2,
     :          0,    1,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=131,140 ) /
     :          0,    1,    0,   -2,    1,
     :         -1,    0,    2,   -2,    2,
     :          0,    0,    0,   -1,    1,
     :         -1,    1,    0,    0,    1,
     :          1,    0,    2,   -1,    2,
     :          1,   -1,    0,    2,    0,
     :          0,    0,    0,    4,    0,
     :          1,    0,    2,    1,    2,
     :          0,    0,    2,    1,    1,
     :          1,    0,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=141,150 ) /
     :         -1,    0,    2,    4,    1,
     :          1,    0,   -2,    0,    1,
     :          1,    1,    2,   -2,    1,
     :          0,    0,    2,    2,    0,
     :         -1,    0,    2,   -1,    1,
     :         -2,    0,    2,    2,    1,
     :          4,    0,    2,    0,    2,
     :          2,   -1,    0,    0,    0,
     :          2,    1,    2,   -2,    2,
     :          0,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=151,160 ) /
     :          1,    0,    4,   -2,    2,
     :         -1,   -1,    0,    0,    1,
     :          0,    1,    0,    2,    1,
     :         -2,    0,    2,    4,    1,
     :          2,    0,    2,    0,    0,
     :          1,    0,    0,    1,    0,
     :         -1,    0,    0,    4,    1,
     :         -1,    0,    4,    0,    1,
     :          2,    0,    2,    2,    1,
     :          0,    0,    2,   -3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=161,170 ) /
     :         -1,   -2,    0,    2,    0,
     :          2,    1,    0,    0,    0,
     :          0,    0,    4,    0,    2,
     :          0,    0,    0,    0,    3,
     :          0,    3,    0,    0,    0,
     :          0,    0,    2,   -4,    1,
     :          0,   -1,    0,    2,    1,
     :          0,    0,    0,    4,    1,
     :         -1,   -1,    2,    4,    2,
     :          1,    0,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=171,180 ) /
     :         -2,    2,    0,    2,    0,
     :         -2,   -1,    2,    0,    1,
     :         -2,    0,    0,    2,    2,
     :         -1,   -1,    2,    0,    2,
     :          0,    0,    4,   -2,    1,
     :          3,    0,    2,   -2,    1,
     :         -2,   -1,    0,    2,    1,
     :          1,    0,    0,   -1,    1,
     :          0,   -2,    0,    2,    0,
     :         -2,    0,    0,    4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=181,190 ) /
     :         -3,    0,    0,    0,    1,
     :          1,    1,    2,    2,    2,
     :          0,    0,    2,    4,    1,
     :          3,    0,    2,    2,    2,
     :         -1,    1,    2,   -2,    1,
     :          2,    0,    0,   -4,    1,
     :          0,    0,    0,   -2,    2,
     :          2,    0,    2,   -4,    1,
     :         -1,    1,    0,    2,    1,
     :          0,    0,    2,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=191,200 ) /
     :          0,   -2,    2,    2,    2,
     :          2,    0,    0,    2,    1,
     :          4,    0,    2,   -2,    2,
     :          2,    0,    0,   -2,    2,
     :          0,    2,    0,    0,    1,
     :          1,    0,    0,   -4,    1,
     :          0,    2,    2,   -2,    1,
     :         -3,    0,    0,    4,    0,
     :         -1,    1,    2,    0,    1,
     :         -1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=201,210 ) /
     :         -1,   -2,    2,    2,    2,
     :         -2,   -1,    2,    4,    2,
     :          1,   -1,    2,    2,    1,
     :         -2,    1,    0,    2,    0,
     :         -2,    1,    2,    0,    1,
     :          2,    1,    0,   -2,    1,
     :         -3,    0,    2,    0,    1,
     :         -2,    0,    2,   -2,    1,
     :         -1,    1,    0,    2,    2,
     :          0,   -1,    2,   -1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=211,220 ) /
     :         -1,    0,    4,   -2,    2,
     :          0,   -2,    2,    0,    2,
     :         -1,    0,    2,    1,    2,
     :          2,    0,    0,    0,    2,
     :          0,    0,    2,    0,    3,
     :         -2,    0,    4,    0,    2,
     :         -1,    0,   -2,    0,    1,
     :         -1,    1,    2,    2,    1,
     :          3,    0,    0,    0,    1,
     :         -1,    0,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=221,230 ) /
     :          2,   -1,    2,    0,    1,
     :          0,    1,    2,    2,    1,
     :          0,   -1,    2,    4,    2,
     :          2,   -1,    2,    2,    2,
     :          0,    2,   -2,    2,    0,
     :         -1,   -1,    2,   -1,    1,
     :          0,   -2,    0,    0,    1,
     :          1,    0,    2,   -4,    2,
     :          1,   -1,    0,   -2,    1,
     :         -1,   -1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=231,240 ) /
     :          1,   -1,    2,   -2,    2,
     :         -2,   -1,    0,    4,    0,
     :         -1,    0,    0,    3,    0,
     :         -2,   -1,    2,    2,    2,
     :          0,    2,    2,    0,    2,
     :          1,    1,    0,    2,    0,
     :          2,    0,    2,   -1,    2,
     :          1,    0,    2,    1,    1,
     :          4,    0,    0,    0,    0,
     :          2,    1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=241,250 ) /
     :          3,   -1,    2,    0,    2,
     :         -2,    2,    0,    2,    1,
     :          1,    0,    2,   -3,    1,
     :          1,    1,    2,   -4,    1,
     :         -1,   -1,    2,   -2,    1,
     :          0,   -1,    0,   -1,    1,
     :          0,   -1,    0,   -2,    1,
     :         -2,    0,    0,    0,    2,
     :         -2,    0,   -2,    2,    0,
     :         -1,    0,   -2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=251,260 ) /
     :          1,   -2,    0,    0,    0,
     :          0,    1,    0,    1,    1,
     :         -1,    2,    0,    2,    0,
     :          1,   -1,    2,   -2,    1,
     :          1,    2,    2,   -2,    2,
     :          2,   -1,    2,   -2,    2,
     :          1,    0,    2,   -1,    1,
     :          2,    1,    2,   -2,    1,
     :         -2,    0,    0,   -2,    1,
     :          1,   -2,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=261,270 ) /
     :          0,    1,    2,    1,    1,
     :          1,    0,    4,   -2,    1,
     :         -2,    0,    4,    2,    2,
     :          1,    1,    2,    1,    2,
     :          1,    0,    0,    4,    0,
     :          1,    0,    2,    2,    0,
     :          2,    0,    2,    1,    2,
     :          3,    1,    2,    0,    2,
     :          4,    0,    2,    0,    1,
     :         -2,   -1,    2,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=271,280 ) /
     :          0,    1,   -2,    2,    1,
     :          1,    0,   -2,    1,    0,
     :          0,   -1,   -2,    2,    1,
     :          2,   -1,    0,   -2,    1,
     :         -1,    0,    2,   -1,    2,
     :          1,    0,    2,   -3,    2,
     :          0,    1,    2,   -2,    3,
     :          0,    0,    2,   -3,    1,
     :         -1,    0,   -2,    2,    1,
     :          0,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=281,290 ) /
     :         -2,    1,    0,    0,    1,
     :         -1,    0,    0,   -1,    1,
     :          2,    0,    2,   -4,    2,
     :          0,    0,    4,   -4,    4,
     :          0,    0,    4,   -4,    2,
     :         -1,   -2,    0,    2,    1,
     :         -2,    0,    0,    3,    0,
     :          1,    0,   -2,    2,    1,
     :         -3,    0,    2,    2,    2,
     :         -3,    0,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=291,300 ) /
     :         -2,    0,    2,    2,    0,
     :          2,   -1,    0,    0,    1,
     :         -2,    1,    2,    2,    2,
     :          1,    1,    0,    1,    0,
     :          0,    1,    4,   -2,    2,
     :         -1,    1,    0,   -2,    1,
     :          0,    0,    0,   -4,    1,
     :          1,   -1,    0,    2,    1,
     :          1,    1,    0,    2,    1,
     :         -1,    2,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=301,310 ) /
     :          3,    1,    2,   -2,    2,
     :          0,   -1,    0,    4,    0,
     :          2,   -1,    0,    2,    0,
     :          0,    0,    4,    0,    1,
     :          2,    0,    4,   -2,    2,
     :         -1,   -1,    2,    4,    1,
     :          1,    0,    0,    4,    1,
     :          1,   -2,    2,    2,    2,
     :          0,    0,    2,    3,    2,
     :         -1,    1,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=311,320 ) /
     :          3,    0,    0,    2,    0,
     :         -1,    0,    4,    2,    2,
     :          1,    1,    2,    2,    1,
     :         -2,    0,    2,    6,    2,
     :          2,    1,    2,    2,    2,
     :         -1,    0,    2,    6,    2,
     :          1,    0,    2,    4,    1,
     :          2,    0,    2,    4,    2,
     :          1,    1,   -2,    1,    0,
     :         -3,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=321,330 ) /
     :          2,    0,   -2,    0,    2,
     :         -1,    0,    0,    1,    2,
     :         -4,    0,    2,    2,    1,
     :         -1,   -1,    0,    1,    0,
     :          0,    0,   -2,    2,    2,
     :          1,    0,    0,   -1,    2,
     :          0,   -1,    2,   -2,    3,
     :         -2,    1,    2,    0,    0,
     :          0,    0,    2,   -2,    4,
     :         -2,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=331,340 ) /
     :         -2,    0,   -2,    4,    0,
     :          0,   -2,   -2,    2,    0,
     :          1,    2,    0,   -2,    1,
     :          3,    0,    0,   -4,    1,
     :         -1,    1,    2,   -2,    2,
     :          1,   -1,    2,   -4,    1,
     :          1,    1,    0,   -2,    2,
     :         -3,    0,    2,    0,    0,
     :         -3,    0,    2,    0,    2,
     :         -2,    0,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=341,350 ) /
     :          0,    0,   -2,    1,    0,
     :         -3,    0,    0,    2,    1,
     :         -1,   -1,   -2,    2,    0,
     :          0,    1,    2,   -4,    1,
     :          2,    1,    0,   -4,    1,
     :          0,    2,    0,   -2,    1,
     :          1,    0,    0,   -3,    1,
     :         -2,    0,    2,   -2,    2,
     :         -2,   -1,    0,    0,    1,
     :         -4,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=351,360 ) /
     :          1,    1,    0,   -4,    1,
     :         -1,    0,    2,   -4,    1,
     :          0,    0,    4,   -4,    1,
     :          0,    3,    2,   -2,    2,
     :         -3,   -1,    0,    4,    0,
     :         -3,    0,    0,    4,    1,
     :          1,   -1,   -2,    2,    0,
     :         -1,   -1,    0,    2,    2,
     :          1,   -2,    0,    0,    1,
     :          1,   -1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=361,370 ) /
     :          0,    0,    0,    1,    2,
     :         -1,   -1,    2,    0,    0,
     :          1,   -2,    2,   -2,    2,
     :          0,   -1,    2,   -1,    1,
     :         -1,    0,    2,    0,    3,
     :          1,    1,    0,    0,    2,
     :         -1,    1,    2,    0,    0,
     :          1,    2,    0,    0,    0,
     :         -1,    2,    2,    0,    2,
     :         -1,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=371,380 ) /
     :          3,    0,    2,   -4,    2,
     :          1,    2,    2,   -2,    1,
     :          1,    0,    4,   -4,    2,
     :         -2,   -1,    0,    4,    1,
     :          0,   -1,    0,    2,    2,
     :         -2,    1,    0,    4,    0,
     :         -2,   -1,    2,    2,    1,
     :          2,    0,   -2,    2,    0,
     :          1,    0,    0,    1,    1,
     :          0,    1,    0,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=381,390 ) /
     :          1,   -1,    2,   -1,    2,
     :         -2,    0,    4,    0,    1,
     :          2,    1,    0,    0,    1,
     :          0,    1,    2,    0,    0,
     :          0,   -1,    4,   -2,    2,
     :          0,    0,    4,   -2,    4,
     :          0,    2,    2,    0,    1,
     :         -3,    0,    0,    6,    0,
     :         -1,   -1,    0,    4,    1,
     :          1,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=391,400 ) /
     :         -1,    0,    0,    4,    2,
     :         -1,   -2,    2,    2,    1,
     :         -1,    0,    0,   -2,    2,
     :          1,    0,   -2,   -2,    1,
     :          0,    0,   -2,   -2,    1,
     :         -2,    0,   -2,    0,    1,
     :          0,    0,    0,    3,    1,
     :          0,    0,    0,    3,    0,
     :         -1,    1,    0,    4,    0,
     :         -1,   -1,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=401,410 ) /
     :         -2,    0,    2,    3,    2,
     :          1,    0,    0,    2,    2,
     :          0,   -1,    2,    1,    2,
     :          3,   -1,    0,    0,    0,
     :          2,    0,    0,    1,    0,
     :          1,   -1,    2,    0,    0,
     :          0,    0,    2,    1,    0,
     :          1,    0,    2,    0,    3,
     :          3,    1,    0,    0,    0,
     :          3,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=411,420 ) /
     :          2,    0,    2,   -1,    1,
     :          1,    1,    2,    0,    0,
     :          0,    0,    4,   -1,    2,
     :          1,    2,    2,    0,    2,
     :         -2,    0,    0,    6,    0,
     :          0,   -1,    0,    4,    1,
     :         -2,   -1,    2,    4,    1,
     :          0,   -2,    2,    2,    1,
     :          0,   -1,    2,    2,    0,
     :         -1,    0,    2,    3,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=421,430 ) /
     :         -2,    1,    2,    4,    2,
     :          2,    0,    0,    2,    2,
     :          2,   -2,    2,    0,    2,
     :         -1,    1,    2,    3,    2,
     :          3,    0,    2,   -1,    2,
     :          4,    0,    2,   -2,    1,
     :         -1,    0,    0,    6,    0,
     :         -1,   -2,    2,    4,    2,
     :         -3,    0,    2,    6,    2,
     :         -1,    0,    2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=431,440 ) /
     :          3,    0,    0,    2,    1,
     :          3,   -1,    2,    0,    1,
     :          3,    0,    2,    0,    0,
     :          1,    0,    4,    0,    2,
     :          5,    0,    2,   -2,    2,
     :          0,   -1,    2,    4,    1,
     :          2,   -1,    2,    2,    1,
     :          0,    1,    2,    4,    2,
     :          1,   -1,    2,    4,    2,
     :          3,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=441,450 ) /
     :          3,    0,    2,    2,    1,
     :          5,    0,    2,    0,    2,
     :          0,    0,    2,    6,    2,
     :          4,    0,    2,    2,    2,
     :          0,   -1,    1,   -1,    1,
     :         -1,    0,    1,    0,    3,
     :          0,   -2,    2,   -2,    3,
     :          1,    0,   -1,    0,    1,
     :          2,   -2,    0,   -2,    1,
     :         -1,    0,    1,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=451,460 ) /
     :         -1,    0,    1,    0,    1,
     :         -1,   -1,    2,   -1,    2,
     :         -2,    2,    0,    2,    2,
     :         -1,    0,    1,    0,    0,
     :         -4,    1,    2,    2,    2,
     :         -3,    0,    2,    1,    1,
     :         -2,   -1,    2,    0,    2,
     :          1,    0,   -2,    1,    1,
     :          2,   -1,   -2,    0,    1,
     :         -4,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=461,470 ) /
     :         -3,    1,    0,    3,    0,
     :         -1,    0,   -1,    2,    0,
     :          0,   -2,    0,    0,    2,
     :          0,   -2,    0,    0,    2,
     :         -3,    0,    0,    3,    0,
     :         -2,   -1,    0,    2,    2,
     :         -1,    0,   -2,    3,    0,
     :         -4,    0,    0,    4,    0,
     :          2,    1,   -2,    0,    1,
     :          2,   -1,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=471,480 ) /
     :          0,    0,    1,   -1,    0,
     :         -1,    2,    0,    1,    0,
     :         -2,    1,    2,    0,    2,
     :          1,    1,    0,   -1,    1,
     :          1,    0,    1,   -2,    1,
     :          0,    2,    0,    0,    2,
     :          1,   -1,    2,   -3,    1,
     :         -1,    1,    2,   -1,    1,
     :         -2,    0,    4,   -2,    2,
     :         -2,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=481,490 ) /
     :         -2,   -2,    0,    2,    1,
     :         -2,    0,   -2,    4,    0,
     :          1,    2,    2,   -4,    1,
     :          1,    1,    2,   -4,    2,
     :         -1,    2,    2,   -2,    1,
     :          2,    0,    0,   -3,    1,
     :         -1,    2,    0,    0,    1,
     :          0,    0,    0,   -2,    0,
     :         -1,   -1,    2,   -2,    2,
     :         -1,    1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=491,500 ) /
     :          0,    0,    0,   -1,    2,
     :         -2,    1,    0,    1,    0,
     :          1,   -2,    0,   -2,    1,
     :          1,    0,   -2,    0,    2,
     :         -3,    1,    0,    2,    0,
     :         -1,    1,   -2,    2,    0,
     :         -1,   -1,    0,    0,    2,
     :         -3,    0,    0,    2,    0,
     :         -3,   -1,    0,    2,    0,
     :          2,    0,    2,   -6,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=501,510 ) /
     :          0,    1,    2,   -4,    2,
     :          2,    0,    0,   -4,    2,
     :         -2,    1,    2,   -2,    1,
     :          0,   -1,    2,   -4,    1,
     :          0,    1,    0,   -2,    2,
     :         -1,    0,    0,   -2,    0,
     :          2,    0,   -2,   -2,    1,
     :         -4,    0,    2,    0,    1,
     :         -1,   -1,    0,   -1,    1,
     :          0,    0,   -2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=511,520 ) /
     :         -3,    0,    0,    1,    0,
     :         -1,    0,   -2,    1,    0,
     :         -2,    0,   -2,    2,    1,
     :          0,    0,   -4,    2,    0,
     :         -2,   -1,   -2,    2,    0,
     :          1,    0,    2,   -6,    1,
     :         -1,    0,    2,   -4,    2,
     :          1,    0,    0,   -4,    2,
     :          2,    1,    2,   -4,    2,
     :          2,    1,    2,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=521,530 ) /
     :          0,    1,    4,   -4,    4,
     :          0,    1,    4,   -4,    2,
     :         -1,   -1,   -2,    4,    0,
     :         -1,   -3,    0,    2,    0,
     :         -1,    0,   -2,    4,    1,
     :         -2,   -1,    0,    3,    0,
     :          0,    0,   -2,    3,    0,
     :         -2,    0,    0,    3,    1,
     :          0,   -1,    0,    1,    0,
     :         -3,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=531,540 ) /
     :          1,    1,   -2,    2,    0,
     :         -1,    1,    0,    2,    2,
     :          1,   -2,    2,   -2,    1,
     :          0,    0,    1,    0,    2,
     :          0,    0,    1,    0,    1,
     :          0,    0,    1,    0,    0,
     :         -1,    2,    0,    2,    1,
     :          0,    0,    2,    0,    2,
     :         -2,    0,    2,    0,    2,
     :          2,    0,    0,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=541,550 ) /
     :          3,    0,    0,   -2,    1,
     :          1,    0,    2,   -2,    3,
     :          1,    2,    0,    0,    1,
     :          2,    0,    2,   -3,    2,
     :         -1,    1,    4,   -2,    2,
     :         -2,   -2,    0,    4,    0,
     :          0,   -3,    0,    2,    0,
     :          0,    0,   -2,    4,    0,
     :         -1,   -1,    0,    3,    0,
     :         -2,    0,    0,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=551,560 ) /
     :         -1,    0,    0,    3,    1,
     :          2,   -2,    0,    0,    0,
     :          1,   -1,    0,    1,    0,
     :         -1,    0,    0,    2,    0,
     :          0,   -2,    2,    0,    1,
     :         -1,    0,    1,    2,    1,
     :         -1,    1,    0,    3,    0,
     :         -1,   -1,    2,    1,    2,
     :          0,   -1,    2,    0,    0,
     :         -2,    1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=561,570 ) /
     :          2,   -2,    2,   -2,    2,
     :          1,    1,    0,    1,    1,
     :          1,    0,    1,    0,    1,
     :          1,    0,    1,    0,    0,
     :          0,    2,    0,    2,    0,
     :          2,   -1,    2,   -2,    1,
     :          0,   -1,    4,   -2,    1,
     :          0,    0,    4,   -2,    3,
     :          0,    1,    4,   -2,    1,
     :          4,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=571,580 ) /
     :          2,    2,    2,   -2,    2,
     :          2,    0,    4,   -4,    2,
     :         -1,   -2,    0,    4,    0,
     :         -1,   -3,    2,    2,    2,
     :         -3,    0,    2,    4,    2,
     :         -3,    0,    2,   -2,    1,
     :         -1,   -1,    0,   -2,    1,
     :         -3,    0,    0,    0,    2,
     :         -3,    0,   -2,    2,    0,
     :          0,    1,    0,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=581,590 ) /
     :         -2,    1,    0,   -2,    1,
     :         -4,    0,    0,    0,    1,
     :         -1,    0,    0,   -4,    1,
     :         -3,    0,    0,   -2,    1,
     :          0,    0,    0,    3,    2,
     :         -1,    1,    0,    4,    1,
     :          1,   -2,    2,    0,    1,
     :          0,    1,    0,    3,    0,
     :         -1,    0,    2,    2,    3,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=591,600 ) /
     :         -2,    0,    2,    2,    2,
     :         -1,    1,    2,    2,    0,
     :          3,    0,    0,    0,    2,
     :          2,    1,    0,    1,    0,
     :          2,   -1,    2,   -1,    2,
     :          0,    0,    2,    0,    1,
     :          0,    0,    3,    0,    3,
     :          0,    0,    3,    0,    2,
     :         -1,    2,    2,    2,    1,
     :         -1,    0,    4,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=601,610 ) /
     :          1,    2,    2,    0,    1,
     :          3,    1,    2,   -2,    1,
     :          1,    1,    4,   -2,    2,
     :         -2,   -1,    0,    6,    0,
     :          0,   -2,    0,    4,    0,
     :         -2,    0,    0,    6,    1,
     :         -2,   -2,    2,    4,    2,
     :          0,   -3,    2,    2,    2,
     :          0,    0,    0,    4,    2,
     :         -1,   -1,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=611,620 ) /
     :         -2,    0,    2,    4,    0,
     :          2,   -1,    0,    2,    1,
     :          1,    0,    0,    3,    0,
     :          0,    1,    0,    4,    1,
     :          0,    1,    0,    4,    0,
     :          1,   -1,    2,    1,    2,
     :          0,    0,    2,    2,    3,
     :          1,    0,    2,    2,    2,
     :         -1,    0,    2,    2,    2,
     :         -2,    0,    4,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=621,630 ) /
     :          2,    1,    0,    2,    1,
     :          2,    1,    0,    2,    0,
     :          2,   -1,    2,    0,    0,
     :          1,    0,    2,    1,    0,
     :          0,    1,    2,    2,    0,
     :          2,    0,    2,    0,    3,
     :          3,    0,    2,    0,    2,
     :          1,    0,    2,    0,    2,
     :          1,    0,    3,    0,    3,
     :          1,    1,    2,    1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=631,640 ) /
     :          0,    2,    2,    2,    2,
     :          2,    1,    2,    0,    0,
     :          2,    0,    4,   -2,    1,
     :          4,    1,    2,   -2,    2,
     :         -1,   -1,    0,    6,    0,
     :         -3,   -1,    2,    6,    2,
     :         -1,    0,    0,    6,    1,
     :         -3,    0,    2,    6,    1,
     :          1,   -1,    0,    4,    1,
     :          1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=641,650 ) /
     :         -2,    0,    2,    5,    2,
     :          1,   -2,    2,    2,    1,
     :          3,   -1,    0,    2,    0,
     :          1,   -1,    2,    2,    0,
     :          0,    0,    2,    3,    1,
     :         -1,    1,    2,    4,    1,
     :          0,    1,    2,    3,    2,
     :         -1,    0,    4,    2,    1,
     :          2,    0,    2,    1,    1,
     :          5,    0,    0,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=651,660 ) /
     :          2,    1,    2,    1,    2,
     :          1,    0,    4,    0,    1,
     :          3,    1,    2,    0,    1,
     :          3,    0,    4,   -2,    2,
     :         -2,   -1,    2,    6,    2,
     :          0,    0,    0,    6,    0,
     :          0,   -2,    2,    4,    2,
     :         -2,    0,    2,    6,    1,
     :          2,    0,    0,    4,    1,
     :          2,    0,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=661,670 ) /
     :          2,   -2,    2,    2,    2,
     :          0,    0,    2,    4,    0,
     :          1,    0,    2,    3,    2,
     :          4,    0,    0,    2,    0,
     :          2,    0,    2,    2,    0,
     :          0,    0,    4,    2,    2,
     :          4,   -1,    2,    0,    2,
     :          3,    0,    2,    1,    2,
     :          2,    1,    2,    2,    1,
     :          4,    1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=671,678 ) /
     :         -1,   -1,    2,    6,    2,
     :         -1,    0,    2,    6,    1,
     :          1,   -1,    2,    4,    1,
     :          1,    1,    2,    4,    2,
     :          3,    1,    2,    2,    2,
     :          5,    0,    2,    0,    1,
     :          2,   -1,    2,    4,    2,
     :          2,    0,    2,    4,    1 /
 
*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*
 
      DATA ( ( CLS(I,J), I=1,6 ), J=  1, 10 ) /
     : -172064161D0, -174666D0,  33386D0, 92052331D0,  9086D0, 15377D0,
     :  -13170906D0,   -1675D0, -13696D0,  5730336D0, -3015D0, -4587D0,
     :   -2276413D0,    -234D0,   2796D0,   978459D0,  -485D0,  1374D0,
     :    2074554D0,     207D0,   -698D0,  -897492D0,   470D0,  -291D0,
     :    1475877D0,   -3633D0,  11817D0,    73871D0,  -184D0, -1924D0,
     :    -516821D0,    1226D0,   -524D0,   224386D0,  -677D0,  -174D0,
     :     711159D0,      73D0,   -872D0,    -6750D0,     0D0,   358D0,
     :    -387298D0,    -367D0,    380D0,   200728D0,    18D0,   318D0,
     :    -301461D0,     -36D0,    816D0,   129025D0,   -63D0,   367D0,
     :     215829D0,    -494D0,    111D0,   -95929D0,   299D0,   132D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 11, 20 ) /
     :     128227D0,     137D0,    181D0,   -68982D0,    -9D0,    39D0,
     :     123457D0,      11D0,     19D0,   -53311D0,    32D0,    -4D0,
     :     156994D0,      10D0,   -168D0,    -1235D0,     0D0,    82D0,
     :      63110D0,      63D0,     27D0,   -33228D0,     0D0,    -9D0,
     :     -57976D0,     -63D0,   -189D0,    31429D0,     0D0,   -75D0,
     :     -59641D0,     -11D0,    149D0,    25543D0,   -11D0,    66D0,
     :     -51613D0,     -42D0,    129D0,    26366D0,     0D0,    78D0,
     :      45893D0,      50D0,     31D0,   -24236D0,   -10D0,    20D0,
     :      63384D0,      11D0,   -150D0,    -1220D0,     0D0,    29D0,
     :     -38571D0,      -1D0,    158D0,    16452D0,   -11D0,    68D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 21, 30 ) /
     :      32481D0,       0D0,      0D0,   -13870D0,     0D0,     0D0,
     :     -47722D0,       0D0,    -18D0,      477D0,     0D0,   -25D0,
     :     -31046D0,      -1D0,    131D0,    13238D0,   -11D0,    59D0,
     :      28593D0,       0D0,     -1D0,   -12338D0,    10D0,    -3D0,
     :      20441D0,      21D0,     10D0,   -10758D0,     0D0,    -3D0,
     :      29243D0,       0D0,    -74D0,     -609D0,     0D0,    13D0,
     :      25887D0,       0D0,    -66D0,     -550D0,     0D0,    11D0,
     :     -14053D0,     -25D0,     79D0,     8551D0,    -2D0,   -45D0,
     :      15164D0,      10D0,     11D0,    -8001D0,     0D0,    -1D0,
     :     -15794D0,      72D0,    -16D0,     6850D0,   -42D0,    -5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 31, 40 ) /
     :      21783D0,       0D0,     13D0,     -167D0,     0D0,    13D0,
     :     -12873D0,     -10D0,    -37D0,     6953D0,     0D0,   -14D0,
     :     -12654D0,      11D0,     63D0,     6415D0,     0D0,    26D0,
     :     -10204D0,       0D0,     25D0,     5222D0,     0D0,    15D0,
     :      16707D0,     -85D0,    -10D0,      168D0,    -1D0,    10D0,
     :      -7691D0,       0D0,     44D0,     3268D0,     0D0,    19D0,
     :     -11024D0,       0D0,    -14D0,      104D0,     0D0,     2D0,
     :       7566D0,     -21D0,    -11D0,    -3250D0,     0D0,    -5D0,
     :      -6637D0,     -11D0,     25D0,     3353D0,     0D0,    14D0,
     :      -7141D0,      21D0,      8D0,     3070D0,     0D0,     4D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 41, 50 ) /
     :      -6302D0,     -11D0,      2D0,     3272D0,     0D0,     4D0,
     :       5800D0,      10D0,      2D0,    -3045D0,     0D0,    -1D0,
     :       6443D0,       0D0,     -7D0,    -2768D0,     0D0,    -4D0,
     :      -5774D0,     -11D0,    -15D0,     3041D0,     0D0,    -5D0,
     :      -5350D0,       0D0,     21D0,     2695D0,     0D0,    12D0,
     :      -4752D0,     -11D0,     -3D0,     2719D0,     0D0,    -3D0,
     :      -4940D0,     -11D0,    -21D0,     2720D0,     0D0,    -9D0,
     :       7350D0,       0D0,     -8D0,      -51D0,     0D0,     4D0,
     :       4065D0,       0D0,      6D0,    -2206D0,     0D0,     1D0,
     :       6579D0,       0D0,    -24D0,     -199D0,     0D0,     2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 51, 60 ) /
     :       3579D0,       0D0,      5D0,    -1900D0,     0D0,     1D0,
     :       4725D0,       0D0,     -6D0,      -41D0,     0D0,     3D0,
     :      -3075D0,       0D0,     -2D0,     1313D0,     0D0,    -1D0,
     :      -2904D0,       0D0,     15D0,     1233D0,     0D0,     7D0,
     :       4348D0,       0D0,    -10D0,      -81D0,     0D0,     2D0,
     :      -2878D0,       0D0,      8D0,     1232D0,     0D0,     4D0,
     :      -4230D0,       0D0,      5D0,      -20D0,     0D0,    -2D0,
     :      -2819D0,       0D0,      7D0,     1207D0,     0D0,     3D0,
     :      -4056D0,       0D0,      5D0,       40D0,     0D0,    -2D0,
     :      -2647D0,       0D0,     11D0,     1129D0,     0D0,     5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 61, 70 ) /
     :      -2294D0,       0D0,    -10D0,     1266D0,     0D0,    -4D0,
     :       2481D0,       0D0,     -7D0,    -1062D0,     0D0,    -3D0,
     :       2179D0,       0D0,     -2D0,    -1129D0,     0D0,    -2D0,
     :       3276D0,       0D0,      1D0,       -9D0,     0D0,     0D0,
     :      -3389D0,       0D0,      5D0,       35D0,     0D0,    -2D0,
     :       3339D0,       0D0,    -13D0,     -107D0,     0D0,     1D0,
     :      -1987D0,       0D0,     -6D0,     1073D0,     0D0,    -2D0,
     :      -1981D0,       0D0,      0D0,      854D0,     0D0,     0D0,
     :       4026D0,       0D0,   -353D0,     -553D0,     0D0,  -139D0,
     :       1660D0,       0D0,     -5D0,     -710D0,     0D0,    -2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 71, 80 ) /
     :      -1521D0,       0D0,      9D0,      647D0,     0D0,     4D0,
     :       1314D0,       0D0,      0D0,     -700D0,     0D0,     0D0,
     :      -1283D0,       0D0,      0D0,      672D0,     0D0,     0D0,
     :      -1331D0,       0D0,      8D0,      663D0,     0D0,     4D0,
     :       1383D0,       0D0,     -2D0,     -594D0,     0D0,    -2D0,
     :       1405D0,       0D0,      4D0,     -610D0,     0D0,     2D0,
     :       1290D0,       0D0,      0D0,     -556D0,     0D0,     0D0,
     :      -1214D0,       0D0,      5D0,      518D0,     0D0,     2D0,
     :       1146D0,       0D0,     -3D0,     -490D0,     0D0,    -1D0,
     :       1019D0,       0D0,     -1D0,     -527D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 81, 90 ) /
     :      -1100D0,       0D0,      9D0,      465D0,     0D0,     4D0,
     :       -970D0,       0D0,      2D0,      496D0,     0D0,     1D0,
     :       1575D0,       0D0,     -6D0,      -50D0,     0D0,     0D0,
     :        934D0,       0D0,     -3D0,     -399D0,     0D0,    -1D0,
     :        922D0,       0D0,     -1D0,     -395D0,     0D0,    -1D0,
     :        815D0,       0D0,     -1D0,     -422D0,     0D0,    -1D0,
     :        834D0,       0D0,      2D0,     -440D0,     0D0,     1D0,
     :       1248D0,       0D0,      0D0,     -170D0,     0D0,     1D0,
     :       1338D0,       0D0,     -5D0,      -39D0,     0D0,     0D0,
     :        716D0,       0D0,     -2D0,     -389D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 91,100 ) /
     :       1282D0,       0D0,     -3D0,      -23D0,     0D0,     1D0,
     :        742D0,       0D0,      1D0,     -391D0,     0D0,     0D0,
     :       1020D0,       0D0,    -25D0,     -495D0,     0D0,   -10D0,
     :        715D0,       0D0,     -4D0,     -326D0,     0D0,     2D0,
     :       -666D0,       0D0,     -3D0,      369D0,     0D0,    -1D0,
     :       -667D0,       0D0,      1D0,      346D0,     0D0,     1D0,
     :       -704D0,       0D0,      0D0,      304D0,     0D0,     0D0,
     :       -694D0,       0D0,      5D0,      294D0,     0D0,     2D0,
     :      -1014D0,       0D0,     -1D0,        4D0,     0D0,    -1D0,
     :       -585D0,       0D0,     -2D0,      316D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=101,110 ) /
     :       -949D0,       0D0,      1D0,        8D0,     0D0,    -1D0,
     :       -595D0,       0D0,      0D0,      258D0,     0D0,     0D0,
     :        528D0,       0D0,      0D0,     -279D0,     0D0,     0D0,
     :       -590D0,       0D0,      4D0,      252D0,     0D0,     2D0,
     :        570D0,       0D0,     -2D0,     -244D0,     0D0,    -1D0,
     :       -502D0,       0D0,      3D0,      250D0,     0D0,     2D0,
     :       -875D0,       0D0,      1D0,       29D0,     0D0,     0D0,
     :       -492D0,       0D0,     -3D0,      275D0,     0D0,    -1D0,
     :        535D0,       0D0,     -2D0,     -228D0,     0D0,    -1D0,
     :       -467D0,       0D0,      1D0,      240D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=111,120 ) /
     :        591D0,       0D0,      0D0,     -253D0,     0D0,     0D0,
     :       -453D0,       0D0,     -1D0,      244D0,     0D0,    -1D0,
     :        766D0,       0D0,      1D0,        9D0,     0D0,     0D0,
     :       -446D0,       0D0,      2D0,      225D0,     0D0,     1D0,
     :       -488D0,       0D0,      2D0,      207D0,     0D0,     1D0,
     :       -468D0,       0D0,      0D0,      201D0,     0D0,     0D0,
     :       -421D0,       0D0,      1D0,      216D0,     0D0,     1D0,
     :        463D0,       0D0,      0D0,     -200D0,     0D0,     0D0,
     :       -673D0,       0D0,      2D0,       14D0,     0D0,     0D0,
     :        658D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=121,130 ) /
     :       -438D0,       0D0,      0D0,      188D0,     0D0,     0D0,
     :       -390D0,       0D0,      0D0,      205D0,     0D0,     0D0,
     :        639D0,     -11D0,     -2D0,      -19D0,     0D0,     0D0,
     :        412D0,       0D0,     -2D0,     -176D0,     0D0,    -1D0,
     :       -361D0,       0D0,      0D0,      189D0,     0D0,     0D0,
     :        360D0,       0D0,     -1D0,     -185D0,     0D0,    -1D0,
     :        588D0,       0D0,     -3D0,      -24D0,     0D0,     0D0,
     :       -578D0,       0D0,      1D0,        5D0,     0D0,     0D0,
     :       -396D0,       0D0,      0D0,      171D0,     0D0,     0D0,
     :        565D0,       0D0,     -1D0,       -6D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=131,140 ) /
     :       -335D0,       0D0,     -1D0,      184D0,     0D0,    -1D0,
     :        357D0,       0D0,      1D0,     -154D0,     0D0,     0D0,
     :        321D0,       0D0,      1D0,     -174D0,     0D0,     0D0,
     :       -301D0,       0D0,     -1D0,      162D0,     0D0,     0D0,
     :       -334D0,       0D0,      0D0,      144D0,     0D0,     0D0,
     :        493D0,       0D0,     -2D0,      -15D0,     0D0,     0D0,
     :        494D0,       0D0,     -2D0,      -19D0,     0D0,     0D0,
     :        337D0,       0D0,     -1D0,     -143D0,     0D0,    -1D0,
     :        280D0,       0D0,     -1D0,     -144D0,     0D0,     0D0,
     :        309D0,       0D0,      1D0,     -134D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=141,150 ) /
     :       -263D0,       0D0,      2D0,      131D0,     0D0,     1D0,
     :        253D0,       0D0,      1D0,     -138D0,     0D0,     0D0,
     :        245D0,       0D0,      0D0,     -128D0,     0D0,     0D0,
     :        416D0,       0D0,     -2D0,      -17D0,     0D0,     0D0,
     :       -229D0,       0D0,      0D0,      128D0,     0D0,     0D0,
     :        231D0,       0D0,      0D0,     -120D0,     0D0,     0D0,
     :       -259D0,       0D0,      2D0,      109D0,     0D0,     1D0,
     :        375D0,       0D0,     -1D0,       -8D0,     0D0,     0D0,
     :        252D0,       0D0,      0D0,     -108D0,     0D0,     0D0,
     :       -245D0,       0D0,      1D0,      104D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=151,160 ) /
     :        243D0,       0D0,     -1D0,     -104D0,     0D0,     0D0,
     :        208D0,       0D0,      1D0,     -112D0,     0D0,     0D0,
     :        199D0,       0D0,      0D0,     -102D0,     0D0,     0D0,
     :       -208D0,       0D0,      1D0,      105D0,     0D0,     0D0,
     :        335D0,       0D0,     -2D0,      -14D0,     0D0,     0D0,
     :       -325D0,       0D0,      1D0,        7D0,     0D0,     0D0,
     :       -187D0,       0D0,      0D0,       96D0,     0D0,     0D0,
     :        197D0,       0D0,     -1D0,     -100D0,     0D0,     0D0,
     :       -192D0,       0D0,      2D0,       94D0,     0D0,     1D0,
     :       -188D0,       0D0,      0D0,       83D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=161,170 ) /
     :        276D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -286D0,       0D0,      1D0,        6D0,     0D0,     0D0,
     :        186D0,       0D0,     -1D0,      -79D0,     0D0,     0D0,
     :       -219D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :        276D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :       -153D0,       0D0,     -1D0,       84D0,     0D0,     0D0,
     :       -156D0,       0D0,      0D0,       81D0,     0D0,     0D0,
     :       -154D0,       0D0,      1D0,       78D0,     0D0,     0D0,
     :       -174D0,       0D0,      1D0,       75D0,     0D0,     0D0,
     :       -163D0,       0D0,      2D0,       69D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=171,180 ) /
     :       -228D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         91D0,       0D0,     -4D0,      -54D0,     0D0,    -2D0,
     :        175D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -159D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        141D0,       0D0,      0D0,      -72D0,     0D0,     0D0,
     :        147D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -132D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        159D0,       0D0,    -28D0,      -54D0,     0D0,    11D0,
     :        213D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        123D0,       0D0,      0D0,      -64D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=181,190 ) /
     :       -118D0,       0D0,     -1D0,       66D0,     0D0,     0D0,
     :        144D0,       0D0,     -1D0,      -61D0,     0D0,     0D0,
     :       -121D0,       0D0,      1D0,       60D0,     0D0,     0D0,
     :       -134D0,       0D0,      1D0,       56D0,     0D0,     1D0,
     :       -105D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       56D0,     0D0,     0D0,
     :        120D0,       0D0,      0D0,      -52D0,     0D0,     0D0,
     :        101D0,       0D0,      0D0,      -54D0,     0D0,     0D0,
     :       -113D0,       0D0,      0D0,       59D0,     0D0,     0D0,
     :       -106D0,       0D0,      0D0,       61D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=191,200 ) /
     :       -129D0,       0D0,      1D0,       55D0,     0D0,     0D0,
     :       -114D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :        113D0,       0D0,     -1D0,      -49D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,       51D0,     0D0,     0D0,
     :       -100D0,       0D0,     -1D0,       56D0,     0D0,     0D0,
     :         87D0,       0D0,      0D0,      -47D0,     0D0,     0D0,
     :        161D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         96D0,       0D0,      0D0,      -50D0,     0D0,     0D0,
     :        151D0,       0D0,     -1D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=201,210 ) /
     :       -104D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :       -110D0,       0D0,      0D0,       48D0,     0D0,     0D0,
     :       -100D0,       0D0,      1D0,       50D0,     0D0,     0D0,
     :         92D0,       0D0,     -5D0,       12D0,     0D0,    -2D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :        -78D0,       0D0,      0D0,       41D0,     0D0,     0D0,
     :        -77D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :          2D0,       0D0,      0D0,       54D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,      -40D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=211,220 ) /
     :        -93D0,       0D0,      0D0,       40D0,     0D0,     0D0,
     :        -83D0,       0D0,     10D0,       40D0,     0D0,    -2D0,
     :         83D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -91D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :        128D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -79D0,       0D0,      0D0,       34D0,     0D0,     0D0,
     :        -83D0,       0D0,      0D0,       47D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -44D0,     0D0,     0D0,
     :         83D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :         91D0,       0D0,      0D0,      -39D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=221,230 ) /
     :        -77D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         68D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :         71D0,       0D0,      0D0,      -31D0,     0D0,     0D0,
     :         62D0,       0D0,      0D0,      -34D0,     0D0,     0D0,
     :        -63D0,       0D0,      0D0,       33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=231,240 ) /
     :        -73D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :        115D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -103D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         63D0,       0D0,      0D0,      -28D0,     0D0,     0D0,
     :         74D0,       0D0,      0D0,      -32D0,     0D0,     0D0,
     :       -103D0,       0D0,     -3D0,        3D0,     0D0,    -1D0,
     :        -69D0,       0D0,      0D0,       30D0,     0D0,     0D0,
     :         57D0,       0D0,      0D0,      -29D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         64D0,       0D0,      0D0,      -33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=241,250 ) /
     :        -63D0,       0D0,      0D0,       26D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       23D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -24D0,     0D0,     0D0,
     :        -48D0,       0D0,      0D0,       25D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :         56D0,       0D0,      0D0,      -25D0,     0D0,     0D0,
     :         88D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :        -75D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=251,260 ) /
     :         85D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         49D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :        -74D0,       0D0,     -3D0,       -1D0,     0D0,    -1D0,
     :        -39D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         51D0,       0D0,      0D0,      -22D0,     0D0,     0D0,
     :        -40D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         41D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :        -42D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,       22D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=261,270 ) /
     :        -42D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         39D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :         46D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -53D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         81D0,       0D0,     -1D0,       -4D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         53D0,       0D0,      0D0,      -23D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=271,280 ) /
     :        -33D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -33D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :        -60D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         48D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,      -14D0,     0D0,     0D0,
     :         38D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=281,290 ) /
     :        -29D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :         44D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=291,300 ) /
     :        -60D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         36D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -37D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=301,310 ) /
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         65D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :        -30D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -31D0,       0D0,      0D0,       13D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=311,320 ) /
     :         49D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       18D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=321,330 ) /
     :        -21D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=331,340 ) /
     :         28D0,       0D0,      0D0,        0D0,     0D0,    -2D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=341,350 ) /
     :        -31D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :         20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=351,360 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         19D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        7D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=361,370 ) /
     :         13D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=371,380 ) /
     :        -26D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=381,390 ) /
     :         22D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         25D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=391,400 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=401,410 ) /
     :         15D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -25D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=411,420 ) /
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         21D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -8D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=421,430 ) /
     :         19D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=431,440 ) /
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :        -24D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=441,450 ) /
     :        -24D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :          0D0,       0D0,  -1988D0,        0D0,     0D0, -1679D0,
     :          0D0,       0D0,    -63D0,        0D0,     0D0,   -27D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     4D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          0D0,       0D0,    364D0,        0D0,     0D0,   176D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=451,460 ) /
     :          0D0,       0D0,  -1044D0,        0D0,     0D0,  -891D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    330D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=461,470 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=471,480 ) /
     :         -5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          0D0,       0D0,    -12D0,        0D0,     0D0,   -10D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=481,490 ) /
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=491,500 ) /
     :         -8D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=501,510 ) /
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=511,520 ) /
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=521,530 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=531,540 ) /
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         10D0,       0D0,     13D0,        6D0,     0D0,    -5D0,
     :          0D0,       0D0,     30D0,        0D0,     0D0,    14D0,
     :          0D0,       0D0,   -162D0,        0D0,     0D0,  -138D0,
     :          0D0,       0D0,     75D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=541,550 ) /
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=551,560 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,     -3D0,        3D0,     0D0,     1D0,
     :          0D0,       0D0,     -3D0,        0D0,     0D0,    -2D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=561,570 ) /
     :         -1D0,       0D0,      3D0,        3D0,     0D0,    -1D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    -13D0,        0D0,     0D0,   -11D0,
     :          3D0,       0D0,      6D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=571,580 ) /
     :          8D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=581,590 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=591,600 ) /
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          0D0,       0D0,    -26D0,        0D0,     0D0,   -11D0,
     :          0D0,       0D0,    -10D0,        0D0,     0D0,    -5D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=601,610 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=611,620 ) /
     :         13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=621,630 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,     -5D0,        0D0,     0D0,    -2D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=631,640 ) /
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=641,650 ) /
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=651,660 ) /
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=661,670 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=671,678 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
 
*
*  Planetary argument multipliers:
*              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre
 
      DATA ( ( NAPL(I,J), I=1,14 ), J=  1, 10 ) /
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 11, 20 ) /
     :         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 21, 30 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2,
     :         2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 31, 40 ) /
     :        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 41, 50 ) /
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0,
     :        -2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :        -1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 51, 60 ) /
     :        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0,
     :        -2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,
     :        -1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 61, 70 ) /
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 71, 80 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 81, 90 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 91,100 ) /
     :        -2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=101,110 ) /
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=111,120 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2,
     :         1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=121,130 ) /
     :         0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=131,140 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=141,150 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=151,160 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1,
     :        -1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=161,170 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,
     :         0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=171,180 ) /
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=181,190 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=191,200 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=201,210 ) /
     :         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=211,220 ) /
     :         0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=221,230 ) /
     :         0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=231,240 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=241,250 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=251,260 ) /
     :         0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2,
     :         0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=261,270 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=271,280 ) /
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2,
     :         0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=281,290 ) /
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=291,300 ) /
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=301,310 ) /
     :         0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=311,320 ) /
     :        -2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=321,330 ) /
     :         0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=331,340 ) /
     :         0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=341,350 ) /
     :         0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=351,360 ) /
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=361,370 ) /
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=371,380 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=381,390 ) /
     :         0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=391,400 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=401,410 ) /
     :         0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=411,420 ) /
     :         0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=421,430 ) /
     :         0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=431,440 ) /
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=441,450 ) /
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=451,460 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=461,470 ) /
     :         0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=471,480 ) /
     :         0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=481,490 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=491,500 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=501,510 ) /
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=511,520 ) /
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=521,530 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=531,540 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=541,550 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=551,560 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=561,570 ) /
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=571,580 ) /
     :         0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=581,590 ) /
     :         0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=591,600 ) /
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=601,610 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=611,620 ) /
     :         0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=621,630 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,
     :         1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=631,640 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=641,650 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=651,660 ) /
     :         0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=661,670 ) /
     :         0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=671,680 ) /
     :         0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=681,687 ) /
     :         1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
 
*
*  Planetary nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, cos), obliquity (sin, cos)
*
 
      DATA ( ( ICPL(I,J), I=1,4 ), J=  1, 10 ) /
     :       1440,          0,          0,          0,
     :         56,       -117,        -42,        -40,
     :        125,        -43,          0,        -54,
     :          0,          5,          0,          0,
     :          3,         -7,         -3,          0,
     :          3,          0,          0,         -2,
     :       -114,          0,          0,         61,
     :       -219,         89,          0,          0,
     :         -3,          0,          0,          0,
     :       -462,       1604,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 11, 20 ) /
     :         99,          0,          0,        -53,
     :         -3,          0,          0,          2,
     :          0,          6,          2,          0,
     :          3,          0,          0,          0,
     :        -12,          0,          0,          0,
     :         14,       -218,        117,          8,
     :         31,       -481,       -257,        -17,
     :       -491,        128,          0,          0,
     :      -3084,       5123,       2735,       1647,
     :      -1444,       2409,      -1286,       -771 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 21, 30 ) /
     :         11,        -24,        -11,         -9,
     :         26,         -9,          0,          0,
     :        103,        -60,          0,          0,
     :          0,        -13,         -7,          0,
     :        -26,        -29,        -16,         14,
     :          9,        -27,        -14,         -5,
     :         12,          0,          0,         -6,
     :         -7,          0,          0,          0,
     :          0,         24,          0,          0,
     :        284,          0,          0,       -151 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 31, 40 ) /
     :        226,        101,          0,          0,
     :          0,         -8,         -2,          0,
     :          0,         -6,         -3,          0,
     :          5,          0,          0,         -3,
     :        -41,        175,         76,         17,
     :          0,         15,          6,          0,
     :        425,        212,       -133,        269,
     :       1200,        598,        319,       -641,
     :        235,        334,          0,          0,
     :         11,        -12,         -7,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 41, 50 ) /
     :          5,         -6,          3,          3,
     :         -5,          0,          0,          3,
     :          6,          0,          0,         -3,
     :         15,          0,          0,          0,
     :         13,          0,          0,         -7,
     :         -6,         -9,          0,          0,
     :        266,        -78,          0,          0,
     :       -460,       -435,       -232,        246,
     :          0,         15,          7,          0,
     :         -3,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 51, 60 ) /
     :          0,        131,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :        -17,        -19,        -10,          9,
     :         -9,        -11,          6,         -5,
     :         -6,          0,          0,          3,
     :        -16,          8,          0,          0,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 61, 70 ) /
     :         11,         24,         11,         -5,
     :         -3,         -4,         -2,          1,
     :          3,          0,          0,         -1,
     :          0,         -8,         -4,          0,
     :          0,          3,          0,          0,
     :          0,          5,          0,          0,
     :          0,          3,          2,          0,
     :         -6,          4,          2,          3,
     :         -3,         -5,          0,          0,
     :         -5,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 71, 80 ) /
     :          4,         24,         13,         -2,
     :        -42,         20,          0,          0,
     :        -10,        233,          0,          0,
     :         -3,          0,          0,          1,
     :         78,        -18,          0,          0,
     :          0,          3,          1,          0,
     :          0,         -3,         -1,          0,
     :          0,         -4,         -2,          1,
     :          0,         -8,         -4,         -1,
     :          0,         -5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 81, 90 ) /
     :         -7,          0,          0,          3,
     :        -14,          8,          3,          6,
     :          0,          8,         -4,          0,
     :          0,         19,         10,          0,
     :         45,        -22,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,          0,          0,
     :          0,          3,          0,          0,
     :          3,          5,          3,         -2,
     :         89,        -16,         -9,        -48 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 91,100 ) /
     :          0,          3,          0,          0,
     :         -3,          7,          4,          2,
     :       -349,        -62,          0,          0,
     :        -15,         22,          0,          0,
     :         -3,          0,          0,          0,
     :        -53,          0,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -8,          0,          0,
     :         15,         -7,         -4,         -8,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=101,110 ) /
     :        -21,        -78,          0,          0,
     :         20,        -70,        -37,        -11,
     :          0,          6,          3,          0,
     :          5,          3,          2,         -2,
     :        -17,         -4,         -2,          9,
     :          0,          6,          3,          0,
     :         32,         15,         -8,         17,
     :        174,         84,         45,        -93,
     :         11,         56,          0,          0,
     :        -66,        -12,         -6,         35 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=111,120 ) /
     :         47,          8,          4,        -25,
     :          0,          8,          4,          0,
     :         10,        -22,        -12,         -5,
     :         -3,          0,          0,          2,
     :        -24,         12,          0,          0,
     :          5,         -6,          0,          0,
     :          3,          0,          0,         -2,
     :          4,          3,          1,         -2,
     :          0,         29,         15,          0,
     :         -5,         -4,         -2,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=121,130 ) /
     :          8,         -3,         -1,         -5,
     :          0,         -3,          0,          0,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -5,          0,          0,          3,
     :         46,         66,         35,        -25,
     :        -14,          7,          0,          0,
     :          0,          3,          2,          0,
     :         -5,          0,          0,          0,
     :        -68,        -34,        -18,         36 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=131,140 ) /
     :          0,         14,          7,          0,
     :         10,         -6,         -3,         -5,
     :         -5,         -4,         -2,          3,
     :         -3,          5,          2,          1,
     :         76,         17,          9,        -41,
     :         84,        298,        159,        -45,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          2,
     :         -3,          0,          0,          1,
     :        -82,        292,        156,         44 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=141,150 ) /
     :        -73,         17,          9,         39,
     :         -9,        -16,          0,          0,
     :          3,          0,         -1,         -2,
     :         -3,          0,          0,          0,
     :         -9,         -5,         -3,          5,
     :       -439,          0,          0,          0,
     :         57,        -28,        -15,        -30,
     :          0,         -6,         -3,          0,
     :         -4,          0,          0,          2,
     :        -40,         57,         30,         21 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=151,160 ) /
     :         23,          7,          3,        -13,
     :        273,         80,         43,       -146,
     :       -449,        430,          0,          0,
     :         -8,        -47,        -25,          4,
     :          6,         47,         25,         -3,
     :          0,         23,         13,          0,
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -2,
     :        -48,       -110,        -59,         26,
     :         51,        114,         61,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=161,170 ) /
     :       -133,          0,          0,         57,
     :          0,          4,          0,          0,
     :        -21,         -6,         -3,         11,
     :          0,         -3,         -1,          0,
     :        -11,        -21,        -11,          6,
     :        -18,       -436,       -233,          9,
     :         35,         -7,          0,          0,
     :          0,          5,          3,          0,
     :         11,         -3,         -1,         -6,
     :         -5,         -3,         -1,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=171,180 ) /
     :        -53,         -9,         -5,         28,
     :          0,          3,          2,          1,
     :          4,          0,          0,         -2,
     :          0,         -4,          0,          0,
     :        -50,        194,        103,         27,
     :        -13,         52,         28,          7,
     :        -91,        248,          0,          0,
     :          6,         49,         26,         -3,
     :         -6,        -47,        -25,          3,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=181,190 ) /
     :         52,         23,         10,        -23,
     :         -3,          0,          0,          1,
     :          0,          5,          3,          0,
     :         -4,          0,          0,          0,
     :         -4,          8,          3,          2,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :          0,          8,          4,          0,
     :          0,          8,          4,          1,
     :         -4,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=191,200 ) /
     :         -4,          0,          0,          0,
     :         -8,          4,          2,          4,
     :          8,         -4,         -2,         -4,
     :          0,         15,          7,          0,
     :       -138,          0,          0,          0,
     :          0,         -7,         -3,          0,
     :          0,         -7,         -3,          0,
     :         54,          0,          0,        -29,
     :          0,         10,          4,          0,
     :         -7,          0,          0,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=201,210 ) /
     :        -37,         35,         19,         20,
     :          0,          4,          0,          0,
     :         -4,          9,          0,          0,
     :          8,          0,          0,         -4,
     :         -9,        -14,         -8,          5,
     :         -3,         -9,         -5,          3,
     :       -145,         47,          0,          0,
     :        -10,         40,         21,          5,
     :         11,        -49,        -26,         -7,
     :      -2150,          0,          0,        932 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=211,220 ) /
     :        -12,          0,          0,          5,
     :         85,          0,          0,        -37,
     :          4,          0,          0,         -2,
     :          3,          0,          0,         -2,
     :        -86,        153,          0,          0,
     :         -6,          9,          5,          3,
     :          9,        -13,         -7,         -5,
     :         -8,         12,          6,          4,
     :        -51,          0,          0,         22,
     :        -11,       -268,       -116,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=221,230 ) /
     :          0,         12,          5,          0,
     :          0,          7,          3,          0,
     :         31,          6,          3,        -17,
     :        140,         27,         14,        -75,
     :         57,         11,          6,        -30,
     :        -14,        -39,          0,          0,
     :          0,         -6,         -2,          0,
     :          4,         15,          8,         -2,
     :          0,          4,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=231,240 ) /
     :          0,         11,          5,          0,
     :          9,          6,          0,          0,
     :         -4,         10,          4,          2,
     :          5,          3,          0,          0,
     :         16,          0,          0,         -9,
     :         -3,          0,          0,          0,
     :          0,          3,          2,         -1,
     :          7,          0,          0,         -3,
     :        -25,         22,          0,          0,
     :         42,        223,        119,        -22 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=241,250 ) /
     :        -27,       -143,        -77,         14,
     :          9,         49,         26,         -5,
     :      -1166,          0,          0,        505,
     :         -5,          0,          0,          2,
     :         -6,          0,          0,          3,
     :         -8,          0,          1,          4,
     :          0,         -4,          0,          0,
     :        117,          0,          0,        -63,
     :         -4,          8,          4,          2,
     :          3,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=251,260 ) /
     :         -5,          0,          0,          2,
     :          0,         31,          0,          0,
     :         -5,          0,          1,          3,
     :          4,          0,          0,         -2,
     :         -4,          0,          0,          2,
     :        -24,        -13,         -6,         10,
     :          3,          0,          0,          0,
     :          0,        -32,        -17,          0,
     :          8,         12,          5,         -3,
     :          3,          0,          0,         -1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=261,270 ) /
     :          7,         13,          0,          0,
     :         -3,         16,          0,          0,
     :         50,          0,          0,        -27,
     :          0,         -5,         -3,          0,
     :         13,          0,          0,          0,
     :          0,          5,          3,          1,
     :         24,          5,          2,        -11,
     :          5,        -11,         -5,         -2,
     :         30,         -3,         -2,        -16,
     :         18,          0,          0,         -9 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=271,280 ) /
     :          8,        614,          0,          0,
     :          3,         -3,         -1,         -2,
     :          6,         17,          9,         -3,
     :         -3,         -9,         -5,          2,
     :          0,          6,          3,         -1,
     :       -127,         21,          9,         55,
     :          3,          5,          0,          0,
     :         -6,        -10,         -4,          3,
     :          5,          0,          0,          0,
     :         16,          9,          4,         -7 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=281,290 ) /
     :          3,          0,          0,         -2,
     :          0,         22,          0,          0,
     :          0,         19,         10,          0,
     :          7,          0,          0,         -4,
     :          0,         -5,         -2,          0,
     :          0,          3,          1,          0,
     :         -9,          3,          1,          4,
     :         17,          0,          0,         -7,
     :          0,         -3,         -2,         -1,
     :        -20,         34,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=291,300 ) /
     :        -10,          0,          1,          5,
     :         -4,          0,          0,          2,
     :         22,        -87,          0,          0,
     :         -4,          0,          0,          2,
     :         -3,         -6,         -2,          1,
     :        -16,         -3,         -1,          7,
     :          0,         -3,         -2,          0,
     :          4,          0,          0,          0,
     :        -68,         39,          0,          0,
     :         27,          0,          0,        -14 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=301,310 ) /
     :          0,         -4,          0,          0,
     :        -25,          0,          0,          0,
     :        -12,         -3,         -2,          6,
     :          3,          0,          0,         -1,
     :          3,         66,         29,         -1,
     :        490,          0,          0,       -213,
     :        -22,         93,         49,         12,
     :         -7,         28,         15,          4,
     :         -3,         13,          7,          2,
     :        -46,         14,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=311,320 ) /
     :         -5,          0,          0,          0,
     :          2,          1,          0,          0,
     :          0,         -3,          0,          0,
     :        -28,          0,          0,         15,
     :          5,          0,          0,         -2,
     :          0,          3,          0,          0,
     :        -11,          0,          0,          5,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          1,
     :         25,        106,         57,        -13 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=321,330 ) /
     :          5,         21,         11,         -3,
     :       1485,          0,          0,          0,
     :         -7,        -32,        -17,          4,
     :          0,          5,          3,          0,
     :         -6,         -3,         -2,          3,
     :         30,         -6,         -2,        -13,
     :         -4,          4,          0,          0,
     :        -19,          0,          0,         10,
     :          0,          4,          2,         -1,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=331,340 ) /
     :          4,          0,          0,         -2,
     :          0,         -3,         -1,          0,
     :         -3,          0,          0,          0,
     :          5,          3,          1,         -2,
     :          0,         11,          0,          0,
     :        118,          0,          0,        -52,
     :          0,         -5,         -3,          0,
     :        -28,         36,          0,          0,
     :          5,         -5,          0,          0,
     :         14,        -59,        -31,         -8 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=341,350 ) /
     :          0,          9,          5,          1,
     :       -458,          0,          0,        198,
     :          0,        -45,        -20,          0,
     :          9,          0,          0,         -5,
     :          0,         -3,          0,          0,
     :          0,         -4,         -2,         -1,
     :         11,          0,          0,         -6,
     :          6,          0,          0,         -2,
     :        -16,         23,          0,          0,
     :          0,         -4,         -2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=351,360 ) /
     :         -5,          0,          0,          2,
     :       -166,        269,          0,          0,
     :         15,          0,          0,         -8,
     :         10,          0,          0,         -4,
     :        -78,         45,          0,          0,
     :          0,         -5,         -2,          0,
     :          7,          0,          0,         -4,
     :         -5,        328,          0,          0,
     :          3,          0,          0,         -2,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=361,370 ) /
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -4,         -2,          0,
     :      -1223,        -26,          0,          0,
     :          0,          7,          3,          0,
     :          3,          0,          0,          0,
     :          0,          3,          2,          0,
     :         -6,         20,          0,          0,
     :       -368,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=371,380 ) /
     :        -75,          0,          0,          0,
     :         11,          0,          0,         -6,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :        -13,        -30,          0,          0,
     :         21,          3,          0,          0,
     :         -3,          0,          0,          1,
     :         -4,          0,          0,          2,
     :          8,        -27,          0,          0,
     :        -19,        -11,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=381,390 ) /
     :         -4,          0,          0,          2,
     :          0,          5,          2,          0,
     :         -6,          0,          0,          2,
     :         -8,          0,          0,          0,
     :         -1,          0,          0,          0,
     :        -14,          0,          0,          6,
     :          6,          0,          0,          0,
     :        -74,          0,          0,         32,
     :          0,         -3,         -1,          0,
     :          4,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=391,400 ) /
     :          8,         11,          0,          0,
     :          0,          3,          2,          0,
     :       -262,          0,          0,        114,
     :          0,         -4,          0,          0,
     :         -7,          0,          0,          4,
     :          0,        -27,        -12,          0,
     :        -19,         -8,         -4,          8,
     :        202,          0,          0,        -87,
     :         -8,         35,         19,          5,
     :          0,          4,          2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=401,410 ) /
     :         16,         -5,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -3,          0,          0,
     :          1,          0,          0,          0,
     :        -35,        -48,        -21,         15,
     :         -3,         -5,         -2,          1,
     :          6,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          0,         -5,          0,          0,
     :         12,         55,         29,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=411,420 ) /
     :          0,          5,          3,          0,
     :       -598,          0,          0,          0,
     :         -3,        -13,         -7,          1,
     :         -5,         -7,         -3,          2,
     :          3,          0,          0,         -1,
     :          5,         -7,          0,          0,
     :          4,          0,          0,         -2,
     :         16,         -6,          0,          0,
     :          8,         -3,          0,          0,
     :          8,        -31,        -16,         -4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=421,430 ) /
     :          0,          3,          1,          0,
     :        113,          0,          0,        -49,
     :          0,        -24,        -10,          0,
     :          4,          0,          0,         -2,
     :         27,          0,          0,          0,
     :         -3,          0,          0,          1,
     :          0,         -4,         -2,          0,
     :          5,          0,          0,         -2,
     :          0,         -3,          0,          0,
     :        -13,          0,          0,          6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=431,440 ) /
     :          5,          0,          0,         -2,
     :        -18,        -10,         -4,          8,
     :         -4,        -28,          0,          0,
     :         -5,          6,          3,          2,
     :         -3,          0,          0,          1,
     :         -5,         -9,         -4,          2,
     :         17,          0,          0,         -7,
     :         11,          4,          0,          0,
     :          0,         -6,         -2,          0,
     :         83,         15,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=441,450 ) /
     :         -4,          0,          0,          2,
     :          0,       -114,        -49,          0,
     :        117,          0,          0,        -51,
     :         -5,         19,         10,          2,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          2,
     :          0,         -3,         -1,          0,
     :          3,          0,          0,          0,
     :          0,         -6,         -2,          0,
     :        393,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=451,460 ) /
     :         -4,         21,         11,          2,
     :         -6,          0,         -1,          3,
     :         -3,          8,          4,          1,
     :          8,          0,          0,          0,
     :         18,        -29,        -13,         -8,
     :          8,         34,         18,         -4,
     :         89,          0,          0,          0,
     :          3,         12,          6,         -1,
     :         54,        -15,         -7,        -24,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=461,470 ) /
     :          3,          0,          0,         -1,
     :          0,         35,          0,          0,
     :       -154,        -30,        -13,         67,
     :         15,          0,          0,          0,
     :          0,          4,          2,          0,
     :          0,          9,          0,          0,
     :         80,        -71,        -31,        -35,
     :          0,        -20,         -9,          0,
     :         11,          5,          2,         -5,
     :         61,        -96,        -42,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=471,480 ) /
     :         14,          9,          4,         -6,
     :        -11,         -6,         -3,          5,
     :          0,         -3,         -1,          0,
     :        123,       -415,       -180,        -53,
     :          0,          0,          0,        -35,
     :         -5,          0,          0,          0,
     :          7,        -32,        -17,         -4,
     :          0,         -9,         -5,          0,
     :          0,         -4,          2,          0,
     :        -89,          0,          0,         38 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=481,490 ) /
     :          0,        -86,        -19,         -6,
     :          0,          0,        -19,          6,
     :       -123,       -416,       -180,         53,
     :          0,         -3,         -1,          0,
     :         12,         -6,         -3,         -5,
     :        -13,          9,          4,          6,
     :          0,        -15,         -7,          0,
     :          3,          0,          0,         -1,
     :        -62,        -97,        -42,         27,
     :        -11,          5,          2,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=491,500 ) /
     :          0,        -19,         -8,          0,
     :         -3,          0,          0,          1,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :        -85,        -70,        -31,         37,
     :        163,        -12,         -5,        -72,
     :        -63,        -16,         -7,         28,
     :        -21,        -32,        -14,          9,
     :          0,         -3,         -1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=501,510 ) /
     :          3,          0,          0,         -2,
     :          0,          8,          0,          0,
     :          3,         10,          4,         -1,
     :          3,          0,          0,         -1,
     :          0,         -7,         -3,          0,
     :          0,         -4,         -2,          0,
     :          6,         19,          0,          0,
     :          5,       -173,        -75,         -2,
     :          0,         -7,         -3,          0,
     :          7,        -12,         -5,         -3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=511,520 ) /
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -1,
     :         74,          0,          0,        -32,
     :         -3,         12,          6,          2,
     :         26,        -14,         -6,        -11,
     :         19,          0,          0,         -8,
     :          6,         24,         13,         -3,
     :         83,          0,          0,          0,
     :          0,        -10,         -5,          0,
     :         11,         -3,         -1,         -5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=521,530 ) /
     :          3,          0,          1,         -1,
     :          3,          0,          0,         -1,
     :         -4,          0,          0,          0,
     :          5,        -23,        -12,         -3,
     :       -339,          0,          0,        147,
     :          0,        -10,         -5,          0,
     :          5,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          0,         -4,         -2,          0,
     :         18,         -3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=531,540 ) /
     :          9,        -11,         -5,         -4,
     :         -8,          0,          0,          4,
     :          3,          0,          0,         -1,
     :          0,          9,          0,          0,
     :          6,         -9,         -4,         -2,
     :         -4,        -12,          0,          0,
     :         67,        -91,        -39,        -29,
     :         30,        -18,         -8,        -13,
     :          0,          0,          0,          0,
     :          0,       -114,        -50,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=541,550 ) /
     :          0,          0,          0,         23,
     :        517,         16,          7,       -224,
     :          0,         -7,         -3,          0,
     :        143,         -3,         -1,        -62,
     :         29,          0,          0,        -13,
     :         -4,          0,          0,          2,
     :         -6,          0,          0,          3,
     :          5,         12,          5,         -2,
     :        -25,          0,          0,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=551,560 ) /
     :          0,          4,          2,          0,
     :        -22,         12,          5,         10,
     :         50,          0,          0,        -22,
     :          0,          7,          4,          0,
     :          0,          3,          1,          0,
     :         -4,          4,          2,          2,
     :         -5,        -11,         -5,          2,
     :          0,          4,          2,          0,
     :          4,         17,          9,         -2,
     :         59,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=561,570 ) /
     :          0,         -4,         -2,          0,
     :         -8,          0,          0,          4,
     :         -3,          0,          0,          0,
     :          4,        -15,         -8,         -2,
     :        370,         -8,          0,       -160,
     :          0,          0,         -3,          0,
     :          0,          3,          1,          0,
     :         -6,          3,          1,          3,
     :          0,          6,          0,          0,
     :        -10,          0,          0,          4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=571,580 ) /
     :          0,          9,          4,          0,
     :          4,         17,          7,         -2,
     :         34,          0,          0,        -15,
     :          0,          5,          3,          0,
     :         -5,          0,          0,          2,
     :        -37,         -7,         -3,         16,
     :          3,         13,          7,         -2,
     :         40,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :       -184,         -3,         -1,         80 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=581,590 ) /
     :         -3,          0,          0,          1,
     :         -3,          0,          0,          0,
     :          0,        -10,         -6,         -1,
     :         31,         -6,          0,        -13,
     :         -3,        -32,        -14,          1,
     :         -7,          0,          0,          3,
     :          0,         -8,         -4,          0,
     :          3,         -4,          0,          0,
     :          0,          4,          0,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=591,600 ) /
     :         19,        -23,        -10,          2,
     :          0,          0,          0,        -10,
     :          0,          3,          2,          0,
     :          0,          9,          5,         -1,
     :         28,          0,          0,          0,
     :          0,         -7,         -4,          0,
     :          8,         -4,          0,         -4,
     :          0,          0,         -2,          0,
     :          0,          3,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=601,610 ) /
     :         -9,          0,          1,          4,
     :          3,         12,          5,         -1,
     :         17,         -3,         -1,          0,
     :          0,          7,          4,          0,
     :         19,          0,          0,          0,
     :          0,         -5,         -3,          0,
     :         14,         -3,          0,         -1,
     :          0,          0,         -1,          0,
     :          0,          0,          0,         -5,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=611,620 ) /
     :         13,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :          2,          9,          4,          3,
     :          0,          0,          0,         -4,
     :          8,          0,          0,          0,
     :          0,          4,          2,          0,
     :          6,          0,          0,         -3,
     :          6,          0,          0,          0,
     :          0,          3,          1,          0,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=621,630 ) /
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          0,
     :          6,          0,          0,          0,
     :          7,          0,          0,          0,
     :         -4,          0,          0,          0,
     :          4,          0,          0,          0,
     :          6,          0,          0,          0,
     :          0,         -4,          0,          0,
     :          0,         -4,          0,          0,
     :          5,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=631,640 ) /
     :         -3,          0,          0,          0,
     :          4,          0,          0,          0,
     :         -5,          0,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :         13,          0,          0,          0,
     :         21,         11,          0,          0,
     :          0,         -5,          0,          0,
     :          0,         -5,         -2,          0,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=641,650 ) /
     :          0,         -5,          0,          0,
     :         -3,          0,          0,          2,
     :         20,         10,          0,          0,
     :        -34,          0,          0,          0,
     :        -19,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :         -6,          0,          0,          3,
     :         -4,          0,          0,          0,
     :          3,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=651,660 ) /
     :          3,          0,          0,          0,
     :          4,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          6,          0,          0,         -3,
     :         -8,          0,          0,          3,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :        126,        -63,        -27,        -55,
     :         -5,          0,          1,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=661,670 ) /
     :         -3,         28,         15,          2,
     :          5,          0,          1,         -2,
     :          0,          9,          4,          1,
     :          0,          9,          4,         -1,
     :       -126,        -63,        -27,         55,
     :          3,          0,          0,         -1,
     :         21,        -11,         -6,        -11,
     :          0,         -4,          0,          0,
     :        -21,        -11,         -6,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=671,680 ) /
     :          0,          3,          1,          0,
     :          8,          0,          0,         -4,
     :         -6,          0,          0,          3,
     :         -3,          0,          0,          1,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          1,
     :         -5,          0,          0,          2,
     :         24,        -12,         -5,        -11,
     :          0,          3,          1,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=681,687 ) /
     :          0,          3,          2,          0,
     :        -24,        -12,         -5,         10,
     :          4,          0,         -1,         -2,
     :         13,          0,          0,         -6,
     :          7,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          3,          0,          0,         -1 /
 
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC
 
*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------
 
*
*  Fundamental (Delaunay) arguments from Simon et al. (1994)
*
 
*  Mean anomaly of the Moon.
      EL  = MOD (         485868.249036D0 +
     :            T*( 1717915923.2178D0 +
     :            T*(         31.8792D0 +
     :            T*(          0.051635D0 +
     :            T*(        - 0.00024470D0 )))), TURNAS ) * DAS2R
 
*  Mean anomaly of the Sun.
      ELP = MOD (        1287104.79305D0 +
     :            T*(  129596581.0481D0 +
     :            T*(        - 0.5532D0 +
     :            T*(          0.000136D0 +
     :            T*(        - 0.00001149D0 )))), TURNAS ) * DAS2R
 
*  Mean argument of the latitude of the Moon.
      F   = MOD (         335779.526232D0 +
     :            T*( 1739527262.8478D0 +
     :            T*(       - 12.7512D0 +
     :            T*(       -  0.001037D0 +
     :            T*(          0.00000417D0 )))), TURNAS ) * DAS2R
 
*  Mean elongation of the Moon from the Sun.
      D   = MOD (        1072260.70369D0 +
     :            T*( 1602961601.2090D0 +
     :            T*(        - 6.3706D0 +
     :            T*(          0.006593D0 +
     :            T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R
 
*  Mean longitude of the ascending node of the Moon.
      OM  = MOD (         450160.398036D0 +
     :            T*(  - 6962890.5431D0 +
     :            T*(          7.4722D0 +
     :            T*(          0.007702D0 +
     :            T*(        - 0.00005939D0 )))), TURNAS ) * DAS2R
 
*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0
 
*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1
 
*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)
 
*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG
 
 100  CONTINUE
 
*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R
 
*  ------------------
*  PLANETARY NUTATION
*  ------------------
 
*  Mean anomaly of the Moon.
      AL   = MOD ( 2.35555598D0 + 8328.6914269554D0 * T, D2PI )
 
*  Mean anomaly of the Sun.
      ALSU = MOD ( 6.24006013D0 + 628.301955D0 * T, D2PI )
 
*  Mean argument of the latitude of the Moon.
      AF   = MOD ( 1.627905234D0 + 8433.466158131D0 * T, D2PI )
 
*  Mean elongation of the Moon from the Sun.
      AD   = MOD ( 5.198466741D0 + 7771.3771468121D0 * T, D2PI )
 
*  Mean longitude of the ascending node of the Moon.
      AOM  = MOD ( 2.18243920D0 - 33.757045D0 * T, D2PI )
 
*  General accumulated precession in longitude.
      APA  = ( 0.02438175D0 + 0.00000538691D0 * T ) * T
 
*  Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
      ALME = MOD ( 4.402608842D0 + 2608.7903141574D0 * T, D2PI )
      ALVE = MOD ( 3.176146697D0 + 1021.3285546211D0 * T, D2PI )
      ALEA = MOD ( 1.753470314D0 +  628.3075849991D0 * T, D2PI )
      ALMA = MOD ( 6.203480913D0 +  334.0612426700D0 * T, D2PI )
      ALJU = MOD ( 0.599546497D0 +   52.9690962641D0 * T, D2PI )
      ALSA = MOD ( 0.874016757D0 +   21.3299104960D0 * T, D2PI )
      ALUR = MOD ( 5.481293871D0 +    7.4781598567D0 * T, D2PI )
      ALNE = MOD ( 5.321159000D0 +    3.8127774000D0 * T, D2PI )
 
*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0
 
*  Summation of planetary nutation series (in reverse order).
      DO 200 I = NPL, 1, -1
 
*     Argument and functions.
         ARG = MOD ( DBLE ( NAPL( 1,I) ) * AL   +
     :               DBLE ( NAPL( 2,I) ) * ALSU +
     :               DBLE ( NAPL( 3,I) ) * AF   +
     :               DBLE ( NAPL( 4,I) ) * AD   +
     :               DBLE ( NAPL( 5,I) ) * AOM  +
     :               DBLE ( NAPL( 6,I) ) * ALME +
     :               DBLE ( NAPL( 7,I) ) * ALVE +
     :               DBLE ( NAPL( 8,I) ) * ALEA +
     :               DBLE ( NAPL( 9,I) ) * ALMA +
     :               DBLE ( NAPL(10,I) ) * ALJU +
     :               DBLE ( NAPL(11,I) ) * ALSA +
     :               DBLE ( NAPL(12,I) ) * ALUR +
     :               DBLE ( NAPL(13,I) ) * ALNE +
     :               DBLE ( NAPL(14,I) ) * APA, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)
 
*     Term.
         DP = DP + DBLE( ICPL(1,I)) * SARG + DBLE( ICPL(2,I)) * CARG
         DE = DE + DBLE( ICPL(3,I)) * SARG + DBLE( ICPL(4,I)) * CARG
 
 200  CONTINUE
 
*  Convert from 0.1 microarcsec units to radians.
      DPSIPL = DP * U2R
      DEPSPL = DE * U2R
 
*  -----
*  TOTAL
*  -----
 
*  Add planetary and luni-solar components.
      DPSI = DPSIPL + DPSILS
      DEPS = DEPSPL + DEPSLS
 
      END
 
 
 
      SUBROUTINE NU2000K ( DATE1, DATE2, DPSI, DEPS )
 
*+
*  - - - - - - - -
*   N U 2 0 0 0 K
*  - - - - - - - -
*
*  Nutation, IAU 2000A model (MHB_2000 without FCN) MODIFIED.  Series
*     truncated for speed of execution, and using Simon et al. (1994)
*     fundamental arguments throughout.  Accuracy, compared to
*     IAU 2000 A series, is 0.1 mas in delta psi and 0.04 mas in
*     delta epsilon and delta psi sin(epsilon) over 6 centuries
*     centered at year 2000 (99% of errors less than these values).
*
*  Modified form of NU2000A, by Pat Wallace, given in subroutine annex
*  to Chapter 5 of IERS Conventions (2003).
*
*  Given:
*     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
*
*  Returned:
*     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
*
*  This revision:  2002 November 25
*                  2004 March 1     (by G. Kaplan)
*
*-----------------------------------------------------------------------
 
      IMPLICIT NONE
 
      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS
 
*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )
 
*  Milliarcseconds to radians
      DOUBLE PRECISION DMAS2R
      PARAMETER ( DMAS2R = DAS2R / 1D3 )
 
*  Arc seconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )
 
*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )
 
*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )
 
*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )
 
*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )
 
*  Miscellaneous
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 DPSILS, DEPSLS,
     :                 ALME, ALVE, ALEA, ALMA, ALJU, ALSA, ALUR, ALNE,
     :                 APA,
     :                 DPSIPL, DEPSPL
      INTEGER I, J
 
*  -------------------------
*  Luni-Solar nutation model
*  -------------------------
 
*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 323 )
 
*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)
 
*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)
 
*  ---------------
*  Planetary terms
*  ---------------
 
*  Number of terms in the planetary nutation model
      INTEGER NPL
      PARAMETER ( NPL = 165 )
 
*  Coefficients for fundamental arguments
      INTEGER NAPL(14,NPL)
 
*  Longitude and obliquity coefficients
      DOUBLE PRECISION CPL(4,NPL)
 
*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------
 
*
*  Luni-Solar argument multipliers:
*               L     L'    F     D     Om
*
      DATA ( ( NALS(I,J), I=1,5 ), J=  1, 20 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2,
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J= 21, 40 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2,
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J= 41, 60 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0,
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J= 61, 80 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2,
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2,
     :         -2,    0,    2,    4,    2,
     :         -1,    0,    4,    0,    2,
     :          2,    0,    2,   -2,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J= 81,100 ) /
     :          2,    0,    2,    2,    2,
     :          1,    0,    0,    2,    1,
     :          3,    0,    0,    0,    0,
     :          3,    0,    2,   -2,    2,
     :          0,    0,    4,   -2,    2,
     :          0,    1,    2,    0,    1,
     :          0,    0,   -2,    2,    1,
     :          0,    0,    2,   -2,    3,
     :         -1,    0,    0,    4,    0,
     :          2,    0,   -2,    0,    1,
     :         -2,    0,    0,    4,    0,
     :         -1,   -1,    0,    2,    1,
     :         -1,    0,    0,    1,    1,
     :          0,    1,    0,    0,    2,
     :          0,    0,   -2,    0,    1,
     :          0,   -1,    2,    0,    1,
     :          0,    0,    2,   -1,    2,
     :          0,    0,    2,    4,    2,
     :         -2,   -1,    0,    2,    0,
     :          1,    1,    0,   -2,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=101,120 ) /
     :         -1,    1,    0,    2,    0,
     :         -1,    1,    0,    1,    2,
     :          1,   -1,    0,    0,    1,
     :          1,   -1,    2,    2,    2,
     :         -1,    1,    2,    2,    2,
     :          3,    0,    2,    0,    1,
     :          0,    1,   -2,    2,    0,
     :         -1,    0,    0,   -2,    1,
     :          0,    1,    2,    2,    2,
     :         -1,   -1,    2,    2,    1,
     :          0,   -1,    0,    0,    2,
     :          1,    0,    2,   -4,    1,
     :         -1,    0,   -2,    2,    0,
     :          0,   -1,    2,    2,    1,
     :          2,   -1,    2,    0,    2,
     :          0,    0,    0,    2,    2,
     :          1,   -1,    2,    0,    1,
     :         -1,    1,    2,    0,    2,
     :          0,    1,    0,    2,    0,
     :          0,   -1,   -2,    2,    0/
      DATA ( ( NALS(I,J), I=1,5 ), J=121,140 ) /
     :          0,    3,    2,   -2,    2,
     :          0,    0,    0,    1,    1,
     :         -1,    0,    2,    2,    0,
     :          2,    1,    2,    0,    2,
     :          1,    1,    0,    0,    1,
     :          1,    1,    2,    0,    1,
     :          2,    0,    0,    2,    0,
     :          1,    0,   -2,    2,    0,
     :         -1,    0,    0,    2,    2,
     :          0,    1,    0,    1,    0,
     :          0,    1,    0,   -2,    1,
     :         -1,    0,    2,   -2,    2,
     :          0,    0,    0,   -1,    1,
     :         -1,    1,    0,    0,    1,
     :          1,    0,    2,   -1,    2,
     :          1,   -1,    0,    2,    0,
     :          0,    0,    0,    4,    0,
     :          1,    0,    2,    1,    2,
     :          0,    0,    2,    1,    1,
     :          1,    0,    0,   -2,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=141,160 ) /
     :         -1,    0,    2,    4,    1,
     :          1,    0,   -2,    0,    1,
     :          1,    1,    2,   -2,    1,
     :          0,    0,    2,    2,    0,
     :         -1,    0,    2,   -1,    1,
     :         -2,    0,    2,    2,    1,
     :          4,    0,    2,    0,    2,
     :          2,   -1,    0,    0,    0,
     :          2,    1,    2,   -2,    2,
     :          0,    1,    2,    1,    2,
     :          1,    0,    4,   -2,    2,
     :         -1,   -1,    0,    0,    1,
     :          0,    1,    0,    2,    1,
     :         -2,    0,    2,    4,    1,
     :          2,    0,    2,    0,    0,
     :          1,    0,    0,    1,    0,
     :         -1,    0,    0,    4,    1,
     :         -1,    0,    4,    0,    1,
     :          2,    0,    2,    2,    1,
     :          0,    0,    2,   -3,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=161,180 ) /
     :         -1,   -2,    0,    2,    0,
     :          2,    1,    0,    0,    0,
     :          0,    0,    4,    0,    2,
     :          0,    0,    0,    0,    3,
     :          0,    3,    0,    0,    0,
     :          0,    0,    2,   -4,    1,
     :          0,   -1,    0,    2,    1,
     :          0,    0,    0,    4,    1,
     :         -1,   -1,    2,    4,    2,
     :          1,    0,    2,    4,    2,
     :         -2,    2,    0,    2,    0,
     :         -2,   -1,    2,    0,    1,
     :         -2,    0,    0,    2,    2,
     :         -1,   -1,    2,    0,    2,
     :          0,    0,    4,   -2,    1,
     :          3,    0,    2,   -2,    1,
     :         -2,   -1,    0,    2,    1,
     :          1,    0,    0,   -1,    1,
     :          0,   -2,    0,    2,    0,
     :         -2,    0,    0,    4,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=181,200 ) /
     :         -3,    0,    0,    0,    1,
     :          1,    1,    2,    2,    2,
     :          0,    0,    2,    4,    1,
     :          3,    0,    2,    2,    2,
     :         -1,    1,    2,   -2,    1,
     :          2,    0,    0,   -4,    1,
     :          0,    0,    0,   -2,    2,
     :          2,    0,    2,   -4,    1,
     :         -1,    1,    0,    2,    1,
     :          0,    0,    2,   -1,    1,
     :          0,   -2,    2,    2,    2,
     :          2,    0,    0,    2,    1,
     :          4,    0,    2,   -2,    2,
     :          2,    0,    0,   -2,    2,
     :          0,    2,    0,    0,    1,
     :          1,    0,    0,   -4,    1,
     :          0,    2,    2,   -2,    1,
     :         -3,    0,    0,    4,    0,
     :         -1,    1,    2,    0,    1,
     :         -1,   -1,    0,    4,    0/
      DATA ( ( NALS(I,J), I=1,5 ), J=201,220 ) /
     :         -1,   -2,    2,    2,    2,
     :         -2,   -1,    2,    4,    2,
     :          1,   -1,    2,    2,    1,
     :         -2,    1,    0,    2,    0,
     :         -2,    1,    2,    0,    1,
     :          2,    1,    0,   -2,    1,
     :         -3,    0,    2,    0,    1,
     :         -2,    0,    2,   -2,    1,
     :         -1,    1,    0,    2,    2,
     :          0,   -1,    2,   -1,    2,
     :         -1,    0,    4,   -2,    2,
     :          0,   -2,    2,    0,    2,
     :         -1,    0,    2,    1,    2,
     :          2,    0,    0,    0,    2,
     :          0,    0,    2,    0,    3,
     :         -2,    0,    4,    0,    2,
     :         -1,    0,   -2,    0,    1,
     :         -1,    1,    2,    2,    1,
     :          3,    0,    0,    0,    1,
     :         -1,    0,    2,    3,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=221,240 ) /
     :          2,   -1,    2,    0,    1,
     :          0,    1,    2,    2,    1,
     :          0,   -1,    2,    4,    2,
     :          2,   -1,    2,    2,    2,
     :          0,    2,   -2,    2,    0,
     :         -1,   -1,    2,   -1,    1,
     :          0,   -2,    0,    0,    1,
     :          1,    0,    2,   -4,    2,
     :          1,   -1,    0,   -2,    1,
     :         -1,   -1,    2,    0,    1,
     :          1,   -1,    2,   -2,    2,
     :         -2,   -1,    0,    4,    0,
     :         -1,    0,    0,    3,    0,
     :         -2,   -1,    2,    2,    2,
     :          0,    2,    2,    0,    2,
     :          1,    1,    0,    2,    0,
     :          2,    0,    2,   -1,    2,
     :          1,    0,    2,    1,    1,
     :          4,    0,    0,    0,    0,
     :          2,    1,    2,    0,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=241,260 ) /
     :          3,   -1,    2,    0,    2,
     :         -2,    2,    0,    2,    1,
     :          1,    0,    2,   -3,    1,
     :          1,    1,    2,   -4,    1,
     :         -1,   -1,    2,   -2,    1,
     :          0,   -1,    0,   -1,    1,
     :          0,   -1,    0,   -2,    1,
     :         -2,    0,    0,    0,    2,
     :         -2,    0,   -2,    2,    0,
     :         -1,    0,   -2,    4,    0,
     :          1,   -2,    0,    0,    0,
     :          0,    1,    0,    1,    1,
     :         -1,    2,    0,    2,    0,
     :          1,   -1,    2,   -2,    1,
     :          1,    2,    2,   -2,    2,
     :          2,   -1,    2,   -2,    2,
     :          1,    0,    2,   -1,    1,
     :          2,    1,    2,   -2,    1,
     :         -2,    0,    0,   -2,    1,
     :          1,   -2,    2,    0,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=261,280 ) /
     :          0,    1,    2,    1,    1,
     :          1,    0,    4,   -2,    1,
     :         -2,    0,    4,    2,    2,
     :          1,    1,    2,    1,    2,
     :          1,    0,    0,    4,    0,
     :          1,    0,    2,    2,    0,
     :          2,    0,    2,    1,    2,
     :          3,    1,    2,    0,    2,
     :          4,    0,    2,    0,    1,
     :         -2,   -1,    2,    0,    0,
     :          0,    1,   -2,    2,    1,
     :          1,    0,   -2,    1,    0,
     :          2,   -1,    0,   -2,    1,
     :         -1,    0,    2,   -1,    2,
     :          1,    0,    2,   -3,    2,
     :          0,    1,    2,   -2,    3,
     :         -1,    0,   -2,    2,    1,
     :          0,    0,    2,   -4,    2,
     :          2,    0,    2,   -4,    2,
     :          0,    0,    4,   -4,    4/
      DATA ( ( NALS(I,J), I=1,5 ), J=281,300 ) /
     :          0,    0,    4,   -4,    2,
     :         -2,    0,    0,    3,    0,
     :          1,    0,   -2,    2,    1,
     :         -3,    0,    2,    2,    2,
     :         -2,    0,    2,    2,    0,
     :          2,   -1,    0,    0,    1,
     :          1,    1,    0,    1,    0,
     :          0,    1,    4,   -2,    2,
     :         -1,    1,    0,   -2,    1,
     :          0,    0,    0,   -4,    1,
     :          1,   -1,    0,    2,    1,
     :          1,    1,    0,    2,    1,
     :         -1,    2,    2,    2,    2,
     :          3,    1,    2,   -2,    2,
     :          0,   -1,    0,    4,    0,
     :          2,   -1,    0,    2,    0,
     :          0,    0,    4,    0,    1,
     :          2,    0,    4,   -2,    2,
     :         -1,   -1,    2,    4,    1,
     :          1,    0,    0,    4,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=301,320 ) /
     :          1,   -2,    2,    2,    2,
     :          0,    0,    2,    3,    2,
     :         -1,    1,    2,    4,    2,
     :          3,    0,    0,    2,    0,
     :         -1,    0,    4,    2,    2,
     :         -2,    0,    2,    6,    2,
     :         -1,    0,    2,    6,    2,
     :          1,    1,   -2,    1,    0,
     :         -1,    0,    0,    1,    2,
     :         -1,   -1,    0,    1,    0,
     :         -2,    0,    0,    1,    0,
     :          0,    0,   -2,    1,    0,
     :          1,   -1,   -2,    2,    0,
     :          1,    2,    0,    0,    0,
     :          3,    0,    2,    0,    0,
     :          0,   -1,    1,   -1,    1,
     :         -1,    0,    1,    0,    3,
     :         -1,    0,    1,    0,    2,
     :         -1,    0,    1,    0,    1,
     :         -1,    0,    1,    0,    0/
      DATA ( ( NALS(I,J), I=1,5 ), J=321,323 ) /
     :          0,    0,    1,    0,    2,
     :          0,    0,    1,    0,    1,
     :          0,    0,    1,    0,    0/
 
*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*
      DATA ( ( CLS(I,J), I=1,6 ), J=  1, 20 ) /
     :-172064161.D0,-174666.D0, 33386.D0, 92052331.D0, 9086.D0,15377.D0,
     : -13170906.D0,  -1675.D0,-13696.D0,  5730336.D0,-3015.D0,-4587.D0,
     :  -2276413.D0,   -234.D0,  2796.D0,   978459.D0, -485.D0, 1374.D0,
     :   2074554.D0,    207.D0,  -698.D0,  -897492.D0,  470.D0, -291.D0,
     :   1475877.D0,  -3633.D0, 11817.D0,    73871.D0, -184.D0,-1924.D0,
     :   -516821.D0,   1226.D0,  -524.D0,   224386.D0, -677.D0, -174.D0,
     :    711159.D0,     73.D0,  -872.D0,    -6750.D0,    0.D0,  358.D0,
     :   -387298.D0,   -367.D0,   380.D0,   200728.D0,   18.D0,  318.D0,
     :   -301461.D0,    -36.D0,   816.D0,   129025.D0,  -63.D0,  367.D0,
     :    215829.D0,   -494.D0,   111.D0,   -95929.D0,  299.D0,  132.D0,
     :    128227.D0,    137.D0,   181.D0,   -68982.D0,   -9.D0,   39.D0,
     :    123457.D0,     11.D0,    19.D0,   -53311.D0,   32.D0,   -4.D0,
     :    156994.D0,     10.D0,  -168.D0,    -1235.D0,    0.D0,   82.D0,
     :     63110.D0,     63.D0,    27.D0,   -33228.D0,    0.D0,   -9.D0,
     :    -57976.D0,    -63.D0,  -189.D0,    31429.D0,    0.D0,  -75.D0,
     :    -59641.D0,    -11.D0,   149.D0,    25543.D0,  -11.D0,   66.D0,
     :    -51613.D0,    -42.D0,   129.D0,    26366.D0,    0.D0,   78.D0,
     :     45893.D0,     50.D0,    31.D0,   -24236.D0,  -10.D0,   20.D0,
     :     63384.D0,     11.D0,  -150.D0,    -1220.D0,    0.D0,   29.D0,
     :    -38571.D0,     -1.D0,   158.D0,    16452.D0,  -11.D0,   68.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 21, 40 ) /
     :     32481.D0,      0.D0,     0.D0,   -13870.D0,    0.D0,    0.D0,
     :    -47722.D0,      0.D0,   -18.D0,      477.D0,    0.D0,  -25.D0,
     :    -31046.D0,     -1.D0,   131.D0,    13238.D0,  -11.D0,   59.D0,
     :     28593.D0,      0.D0,    -1.D0,   -12338.D0,   10.D0,   -3.D0,
     :     20441.D0,     21.D0,    10.D0,   -10758.D0,    0.D0,   -3.D0,
     :     29243.D0,      0.D0,   -74.D0,     -609.D0,    0.D0,   13.D0,
     :     25887.D0,      0.D0,   -66.D0,     -550.D0,    0.D0,   11.D0,
     :    -14053.D0,    -25.D0,    79.D0,     8551.D0,   -2.D0,  -45.D0,
     :     15164.D0,     10.D0,    11.D0,    -8001.D0,    0.D0,   -1.D0,
     :    -15794.D0,     72.D0,   -16.D0,     6850.D0,  -42.D0,   -5.D0,
     :     21783.D0,      0.D0,    13.D0,     -167.D0,    0.D0,   13.D0,
     :    -12873.D0,    -10.D0,   -37.D0,     6953.D0,    0.D0,  -14.D0,
     :    -12654.D0,     11.D0,    63.D0,     6415.D0,    0.D0,   26.D0,
     :    -10204.D0,      0.D0,    25.D0,     5222.D0,    0.D0,   15.D0,
     :     16707.D0,    -85.D0,   -10.D0,      168.D0,   -1.D0,   10.D0,
     :     -7691.D0,      0.D0,    44.D0,     3268.D0,    0.D0,   19.D0,
     :    -11024.D0,      0.D0,   -14.D0,      104.D0,    0.D0,    2.D0,
     :      7566.D0,    -21.D0,   -11.D0,    -3250.D0,    0.D0,   -5.D0,
     :     -6637.D0,    -11.D0,    25.D0,     3353.D0,    0.D0,   14.D0,
     :     -7141.D0,     21.D0,     8.D0,     3070.D0,    0.D0,    4.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 41, 60 ) /
     :     -6302.D0,    -11.D0,     2.D0,     3272.D0,    0.D0,    4.D0,
     :      5800.D0,     10.D0,     2.D0,    -3045.D0,    0.D0,   -1.D0,
     :      6443.D0,      0.D0,    -7.D0,    -2768.D0,    0.D0,   -4.D0,
     :     -5774.D0,    -11.D0,   -15.D0,     3041.D0,    0.D0,   -5.D0,
     :     -5350.D0,      0.D0,    21.D0,     2695.D0,    0.D0,   12.D0,
     :     -4752.D0,    -11.D0,    -3.D0,     2719.D0,    0.D0,   -3.D0,
     :     -4940.D0,    -11.D0,   -21.D0,     2720.D0,    0.D0,   -9.D0,
     :      7350.D0,      0.D0,    -8.D0,      -51.D0,    0.D0,    4.D0,
     :      4065.D0,      0.D0,     6.D0,    -2206.D0,    0.D0,    1.D0,
     :      6579.D0,      0.D0,   -24.D0,     -199.D0,    0.D0,    2.D0,
     :      3579.D0,      0.D0,     5.D0,    -1900.D0,    0.D0,    1.D0,
     :      4725.D0,      0.D0,    -6.D0,      -41.D0,    0.D0,    3.D0,
     :     -3075.D0,      0.D0,    -2.D0,     1313.D0,    0.D0,   -1.D0,
     :     -2904.D0,      0.D0,    15.D0,     1233.D0,    0.D0,    7.D0,
     :      4348.D0,      0.D0,   -10.D0,      -81.D0,    0.D0,    2.D0,
     :     -2878.D0,      0.D0,     8.D0,     1232.D0,    0.D0,    4.D0,
     :     -4230.D0,      0.D0,     5.D0,      -20.D0,    0.D0,   -2.D0,
     :     -2819.D0,      0.D0,     7.D0,     1207.D0,    0.D0,    3.D0,
     :     -4056.D0,      0.D0,     5.D0,       40.D0,    0.D0,   -2.D0,
     :     -2647.D0,      0.D0,    11.D0,     1129.D0,    0.D0,    5.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 61, 80 ) /
     :     -2294.D0,      0.D0,   -10.D0,     1266.D0,    0.D0,   -4.D0,
     :      2481.D0,      0.D0,    -7.D0,    -1062.D0,    0.D0,   -3.D0,
     :      2179.D0,      0.D0,    -2.D0,    -1129.D0,    0.D0,   -2.D0,
     :      3276.D0,      0.D0,     1.D0,       -9.D0,    0.D0,    0.D0,
     :     -3389.D0,      0.D0,     5.D0,       35.D0,    0.D0,   -2.D0,
     :      3339.D0,      0.D0,   -13.D0,     -107.D0,    0.D0,    1.D0,
     :     -1987.D0,      0.D0,    -6.D0,     1073.D0,    0.D0,   -2.D0,
     :     -1981.D0,      0.D0,     0.D0,      854.D0,    0.D0,    0.D0,
     :      4026.D0,      0.D0,  -353.D0,     -553.D0,    0.D0, -139.D0,
     :      1660.D0,      0.D0,    -5.D0,     -710.D0,    0.D0,   -2.D0,
     :     -1521.D0,      0.D0,     9.D0,      647.D0,    0.D0,    4.D0,
     :      1314.D0,      0.D0,     0.D0,     -700.D0,    0.D0,    0.D0,
     :     -1283.D0,      0.D0,     0.D0,      672.D0,    0.D0,    0.D0,
     :     -1331.D0,      0.D0,     8.D0,      663.D0,    0.D0,    4.D0,
     :      1383.D0,      0.D0,    -2.D0,     -594.D0,    0.D0,   -2.D0,
     :      1405.D0,      0.D0,     4.D0,     -610.D0,    0.D0,    2.D0,
     :      1290.D0,      0.D0,     0.D0,     -556.D0,    0.D0,    0.D0,
     :     -1214.D0,      0.D0,     5.D0,      518.D0,    0.D0,    2.D0,
     :      1146.D0,      0.D0,    -3.D0,     -490.D0,    0.D0,   -1.D0,
     :      1019.D0,      0.D0,    -1.D0,     -527.D0,    0.D0,   -1.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 81,100 ) /
     :     -1100.D0,      0.D0,     9.D0,      465.D0,    0.D0,    4.D0,
     :      -970.D0,      0.D0,     2.D0,      496.D0,    0.D0,    1.D0,
     :      1575.D0,      0.D0,    -6.D0,      -50.D0,    0.D0,    0.D0,
     :       934.D0,      0.D0,    -3.D0,     -399.D0,    0.D0,   -1.D0,
     :       922.D0,      0.D0,    -1.D0,     -395.D0,    0.D0,   -1.D0,
     :       815.D0,      0.D0,    -1.D0,     -422.D0,    0.D0,   -1.D0,
     :       834.D0,      0.D0,     2.D0,     -440.D0,    0.D0,    1.D0,
     :      1248.D0,      0.D0,     0.D0,     -170.D0,    0.D0,    1.D0,
     :      1338.D0,      0.D0,    -5.D0,      -39.D0,    0.D0,    0.D0,
     :       716.D0,      0.D0,    -2.D0,     -389.D0,    0.D0,   -1.D0,
     :      1282.D0,      0.D0,    -3.D0,      -23.D0,    0.D0,    1.D0,
     :       742.D0,      0.D0,     1.D0,     -391.D0,    0.D0,    0.D0,
     :      1020.D0,      0.D0,   -25.D0,     -495.D0,    0.D0,  -10.D0,
     :       715.D0,      0.D0,    -4.D0,     -326.D0,    0.D0,    2.D0,
     :      -666.D0,      0.D0,    -3.D0,      369.D0,    0.D0,   -1.D0,
     :      -667.D0,      0.D0,     1.D0,      346.D0,    0.D0,    1.D0,
     :      -704.D0,      0.D0,     0.D0,      304.D0,    0.D0,    0.D0,
     :      -694.D0,      0.D0,     5.D0,      294.D0,    0.D0,    2.D0,
     :     -1014.D0,      0.D0,    -1.D0,        4.D0,    0.D0,   -1.D0,
     :      -585.D0,      0.D0,    -2.D0,      316.D0,    0.D0,   -1.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=101,120 ) /
     :      -949.D0,      0.D0,     1.D0,        8.D0,    0.D0,   -1.D0,
     :      -595.D0,      0.D0,     0.D0,      258.D0,    0.D0,    0.D0,
     :       528.D0,      0.D0,     0.D0,     -279.D0,    0.D0,    0.D0,
     :      -590.D0,      0.D0,     4.D0,      252.D0,    0.D0,    2.D0,
     :       570.D0,      0.D0,    -2.D0,     -244.D0,    0.D0,   -1.D0,
     :      -502.D0,      0.D0,     3.D0,      250.D0,    0.D0,    2.D0,
     :      -875.D0,      0.D0,     1.D0,       29.D0,    0.D0,    0.D0,
     :      -492.D0,      0.D0,    -3.D0,      275.D0,    0.D0,   -1.D0,
     :       535.D0,      0.D0,    -2.D0,     -228.D0,    0.D0,   -1.D0,
     :      -467.D0,      0.D0,     1.D0,      240.D0,    0.D0,    1.D0,
     :       591.D0,      0.D0,     0.D0,     -253.D0,    0.D0,    0.D0,
     :      -453.D0,      0.D0,    -1.D0,      244.D0,    0.D0,   -1.D0,
     :       766.D0,      0.D0,     1.D0,        9.D0,    0.D0,    0.D0,
     :      -446.D0,      0.D0,     2.D0,      225.D0,    0.D0,    1.D0,
     :      -488.D0,      0.D0,     2.D0,      207.D0,    0.D0,    1.D0,
     :      -468.D0,      0.D0,     0.D0,      201.D0,    0.D0,    0.D0,
     :      -421.D0,      0.D0,     1.D0,      216.D0,    0.D0,    1.D0,
     :       463.D0,      0.D0,     0.D0,     -200.D0,    0.D0,    0.D0,
     :      -673.D0,      0.D0,     2.D0,       14.D0,    0.D0,    0.D0,
     :       658.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=121,140 ) /
     :      -438.D0,      0.D0,     0.D0,      188.D0,    0.D0,    0.D0,
     :      -390.D0,      0.D0,     0.D0,      205.D0,    0.D0,    0.D0,
     :       639.D0,    -11.D0,    -2.D0,      -19.D0,    0.D0,    0.D0,
     :       412.D0,      0.D0,    -2.D0,     -176.D0,    0.D0,   -1.D0,
     :      -361.D0,      0.D0,     0.D0,      189.D0,    0.D0,    0.D0,
     :       360.D0,      0.D0,    -1.D0,     -185.D0,    0.D0,   -1.D0,
     :       588.D0,      0.D0,    -3.D0,      -24.D0,    0.D0,    0.D0,
     :      -578.D0,      0.D0,     1.D0,        5.D0,    0.D0,    0.D0,
     :      -396.D0,      0.D0,     0.D0,      171.D0,    0.D0,    0.D0,
     :       565.D0,      0.D0,    -1.D0,       -6.D0,    0.D0,    0.D0,
     :      -335.D0,      0.D0,    -1.D0,      184.D0,    0.D0,   -1.D0,
     :       357.D0,      0.D0,     1.D0,     -154.D0,    0.D0,    0.D0,
     :       321.D0,      0.D0,     1.D0,     -174.D0,    0.D0,    0.D0,
     :      -301.D0,      0.D0,    -1.D0,      162.D0,    0.D0,    0.D0,
     :      -334.D0,      0.D0,     0.D0,      144.D0,    0.D0,    0.D0,
     :       493.D0,      0.D0,    -2.D0,      -15.D0,    0.D0,    0.D0,
     :       494.D0,      0.D0,    -2.D0,      -19.D0,    0.D0,    0.D0,
     :       337.D0,      0.D0,    -1.D0,     -143.D0,    0.D0,   -1.D0,
     :       280.D0,      0.D0,    -1.D0,     -144.D0,    0.D0,    0.D0,
     :       309.D0,      0.D0,     1.D0,     -134.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=141,160 ) /
     :      -263.D0,      0.D0,     2.D0,      131.D0,    0.D0,    1.D0,
     :       253.D0,      0.D0,     1.D0,     -138.D0,    0.D0,    0.D0,
     :       245.D0,      0.D0,     0.D0,     -128.D0,    0.D0,    0.D0,
     :       416.D0,      0.D0,    -2.D0,      -17.D0,    0.D0,    0.D0,
     :      -229.D0,      0.D0,     0.D0,      128.D0,    0.D0,    0.D0,
     :       231.D0,      0.D0,     0.D0,     -120.D0,    0.D0,    0.D0,
     :      -259.D0,      0.D0,     2.D0,      109.D0,    0.D0,    1.D0,
     :       375.D0,      0.D0,    -1.D0,       -8.D0,    0.D0,    0.D0,
     :       252.D0,      0.D0,     0.D0,     -108.D0,    0.D0,    0.D0,
     :      -245.D0,      0.D0,     1.D0,      104.D0,    0.D0,    0.D0,
     :       243.D0,      0.D0,    -1.D0,     -104.D0,    0.D0,    0.D0,
     :       208.D0,      0.D0,     1.D0,     -112.D0,    0.D0,    0.D0,
     :       199.D0,      0.D0,     0.D0,     -102.D0,    0.D0,    0.D0,
     :      -208.D0,      0.D0,     1.D0,      105.D0,    0.D0,    0.D0,
     :       335.D0,      0.D0,    -2.D0,      -14.D0,    0.D0,    0.D0,
     :      -325.D0,      0.D0,     1.D0,        7.D0,    0.D0,    0.D0,
     :      -187.D0,      0.D0,     0.D0,       96.D0,    0.D0,    0.D0,
     :       197.D0,      0.D0,    -1.D0,     -100.D0,    0.D0,    0.D0,
     :      -192.D0,      0.D0,     2.D0,       94.D0,    0.D0,    1.D0,
     :      -188.D0,      0.D0,     0.D0,       83.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=161,180 ) /
     :       276.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :      -286.D0,      0.D0,     1.D0,        6.D0,    0.D0,    0.D0,
     :       186.D0,      0.D0,    -1.D0,      -79.D0,    0.D0,    0.D0,
     :      -219.D0,      0.D0,     0.D0,       43.D0,    0.D0,    0.D0,
     :       276.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :      -153.D0,      0.D0,    -1.D0,       84.D0,    0.D0,    0.D0,
     :      -156.D0,      0.D0,     0.D0,       81.D0,    0.D0,    0.D0,
     :      -154.D0,      0.D0,     1.D0,       78.D0,    0.D0,    0.D0,
     :      -174.D0,      0.D0,     1.D0,       75.D0,    0.D0,    0.D0,
     :      -163.D0,      0.D0,     2.D0,       69.D0,    0.D0,    1.D0,
     :      -228.D0,      0.D0,     0.D0,        1.D0,    0.D0,    0.D0,
     :        91.D0,      0.D0,    -4.D0,      -54.D0,    0.D0,   -2.D0,
     :       175.D0,      0.D0,     0.D0,      -75.D0,    0.D0,    0.D0,
     :      -159.D0,      0.D0,     0.D0,       69.D0,    0.D0,    0.D0,
     :       141.D0,      0.D0,     0.D0,      -72.D0,    0.D0,    0.D0,
     :       147.D0,      0.D0,     0.D0,      -75.D0,    0.D0,    0.D0,
     :      -132.D0,      0.D0,     0.D0,       69.D0,    0.D0,    0.D0,
     :       159.D0,      0.D0,   -28.D0,      -54.D0,    0.D0,   11.D0,
     :       213.D0,      0.D0,     0.D0,       -4.D0,    0.D0,    0.D0,
     :       123.D0,      0.D0,     0.D0,      -64.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=181,200 ) /
     :      -118.D0,      0.D0,    -1.D0,       66.D0,    0.D0,    0.D0,
     :       144.D0,      0.D0,    -1.D0,      -61.D0,    0.D0,    0.D0,
     :      -121.D0,      0.D0,     1.D0,       60.D0,    0.D0,    0.D0,
     :      -134.D0,      0.D0,     1.D0,       56.D0,    0.D0,    1.D0,
     :      -105.D0,      0.D0,     0.D0,       57.D0,    0.D0,    0.D0,
     :      -102.D0,      0.D0,     0.D0,       56.D0,    0.D0,    0.D0,
     :       120.D0,      0.D0,     0.D0,      -52.D0,    0.D0,    0.D0,
     :       101.D0,      0.D0,     0.D0,      -54.D0,    0.D0,    0.D0,
     :      -113.D0,      0.D0,     0.D0,       59.D0,    0.D0,    0.D0,
     :      -106.D0,      0.D0,     0.D0,       61.D0,    0.D0,    0.D0,
     :      -129.D0,      0.D0,     1.D0,       55.D0,    0.D0,    0.D0,
     :      -114.D0,      0.D0,     0.D0,       57.D0,    0.D0,    0.D0,
     :       113.D0,      0.D0,    -1.D0,      -49.D0,    0.D0,    0.D0,
     :      -102.D0,      0.D0,     0.D0,       44.D0,    0.D0,    0.D0,
     :       -94.D0,      0.D0,     0.D0,       51.D0,    0.D0,    0.D0,
     :      -100.D0,      0.D0,    -1.D0,       56.D0,    0.D0,    0.D0,
     :        87.D0,      0.D0,     0.D0,      -47.D0,    0.D0,    0.D0,
     :       161.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :        96.D0,      0.D0,     0.D0,      -50.D0,    0.D0,    0.D0,
     :       151.D0,      0.D0,    -1.D0,       -5.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=201,220 ) /
     :      -104.D0,      0.D0,     0.D0,       44.D0,    0.D0,    0.D0,
     :      -110.D0,      0.D0,     0.D0,       48.D0,    0.D0,    0.D0,
     :      -100.D0,      0.D0,     1.D0,       50.D0,    0.D0,    0.D0,
     :        92.D0,      0.D0,    -5.D0,       12.D0,    0.D0,   -2.D0,
     :        82.D0,      0.D0,     0.D0,      -45.D0,    0.D0,    0.D0,
     :        82.D0,      0.D0,     0.D0,      -45.D0,    0.D0,    0.D0,
     :       -78.D0,      0.D0,     0.D0,       41.D0,    0.D0,    0.D0,
     :       -77.D0,      0.D0,     0.D0,       43.D0,    0.D0,    0.D0,
     :         2.D0,      0.D0,     0.D0,       54.D0,    0.D0,    0.D0,
     :        94.D0,      0.D0,     0.D0,      -40.D0,    0.D0,    0.D0,
     :       -93.D0,      0.D0,     0.D0,       40.D0,    0.D0,    0.D0,
     :       -83.D0,      0.D0,    10.D0,       40.D0,    0.D0,   -2.D0,
     :        83.D0,      0.D0,     0.D0,      -36.D0,    0.D0,    0.D0,
     :       -91.D0,      0.D0,     0.D0,       39.D0,    0.D0,    0.D0,
     :       128.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :       -79.D0,      0.D0,     0.D0,       34.D0,    0.D0,    0.D0,
     :       -83.D0,      0.D0,     0.D0,       47.D0,    0.D0,    0.D0,
     :        84.D0,      0.D0,     0.D0,      -44.D0,    0.D0,    0.D0,
     :        83.D0,      0.D0,     0.D0,      -43.D0,    0.D0,    0.D0,
     :        91.D0,      0.D0,     0.D0,      -39.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=221,240 ) /
     :       -77.D0,      0.D0,     0.D0,       39.D0,    0.D0,    0.D0,
     :        84.D0,      0.D0,     0.D0,      -43.D0,    0.D0,    0.D0,
     :       -92.D0,      0.D0,     1.D0,       39.D0,    0.D0,    0.D0,
     :       -92.D0,      0.D0,     1.D0,       39.D0,    0.D0,    0.D0,
     :       -94.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        68.D0,      0.D0,     0.D0,      -36.D0,    0.D0,    0.D0,
     :       -61.D0,      0.D0,     0.D0,       32.D0,    0.D0,    0.D0,
     :        71.D0,      0.D0,     0.D0,      -31.D0,    0.D0,    0.D0,
     :        62.D0,      0.D0,     0.D0,      -34.D0,    0.D0,    0.D0,
     :       -63.D0,      0.D0,     0.D0,       33.D0,    0.D0,    0.D0,
     :       -73.D0,      0.D0,     0.D0,       32.D0,    0.D0,    0.D0,
     :       115.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :      -103.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :        63.D0,      0.D0,     0.D0,      -28.D0,    0.D0,    0.D0,
     :        74.D0,      0.D0,     0.D0,      -32.D0,    0.D0,    0.D0,
     :      -103.D0,      0.D0,    -3.D0,        3.D0,    0.D0,   -1.D0,
     :       -69.D0,      0.D0,     0.D0,       30.D0,    0.D0,    0.D0,
     :        57.D0,      0.D0,     0.D0,      -29.D0,    0.D0,    0.D0,
     :        94.D0,      0.D0,     0.D0,       -4.D0,    0.D0,    0.D0,
     :        64.D0,      0.D0,     0.D0,      -33.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=241,260 ) /
     :       -63.D0,      0.D0,     0.D0,       26.D0,    0.D0,    0.D0,
     :       -38.D0,      0.D0,     0.D0,       20.D0,    0.D0,    0.D0,
     :       -43.D0,      0.D0,     0.D0,       24.D0,    0.D0,    0.D0,
     :       -45.D0,      0.D0,     0.D0,       23.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,      -24.D0,    0.D0,    0.D0,
     :       -48.D0,      0.D0,     0.D0,       25.D0,    0.D0,    0.D0,
     :        45.D0,      0.D0,     0.D0,      -26.D0,    0.D0,    0.D0,
     :        56.D0,      0.D0,     0.D0,      -25.D0,    0.D0,    0.D0,
     :        88.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :       -75.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        85.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        49.D0,      0.D0,     0.D0,      -26.D0,    0.D0,    0.D0,
     :       -74.D0,      0.D0,    -3.D0,       -1.D0,    0.D0,   -1.D0,
     :       -39.D0,      0.D0,     0.D0,       21.D0,    0.D0,    0.D0,
     :        45.D0,      0.D0,     0.D0,      -20.D0,    0.D0,    0.D0,
     :        51.D0,      0.D0,     0.D0,      -22.D0,    0.D0,    0.D0,
     :       -40.D0,      0.D0,     0.D0,       21.D0,    0.D0,    0.D0,
     :        41.D0,      0.D0,     0.D0,      -21.D0,    0.D0,    0.D0,
     :       -42.D0,      0.D0,     0.D0,       24.D0,    0.D0,    0.D0,
     :       -51.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=261,280 ) /
     :       -42.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0,
     :        39.D0,      0.D0,     0.D0,      -21.D0,    0.D0,    0.D0,
     :        46.D0,      0.D0,     0.D0,      -18.D0,    0.D0,    0.D0,
     :       -53.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0,
     :        82.D0,      0.D0,     0.D0,       -4.D0,    0.D0,    0.D0,
     :        81.D0,      0.D0,    -1.D0,       -4.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,      -19.D0,    0.D0,    0.D0,
     :        53.D0,      0.D0,     0.D0,      -23.D0,    0.D0,    0.D0,
     :       -45.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0,
     :       -44.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :       -33.D0,      0.D0,     0.D0,       16.D0,    0.D0,    0.D0,
     :       -61.D0,      0.D0,     0.D0,        1.D0,    0.D0,    0.D0,
     :       -38.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :       -33.D0,      0.D0,     0.D0,       21.D0,    0.D0,    0.D0,
     :       -60.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        48.D0,      0.D0,     0.D0,      -10.D0,    0.D0,    0.D0,
     :        38.D0,      0.D0,     0.D0,      -20.D0,    0.D0,    0.D0,
     :        31.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :       -32.D0,      0.D0,     0.D0,       15.D0,    0.D0,    0.D0,
     :        45.D0,      0.D0,     0.D0,       -8.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=281,300 ) /
     :       -44.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :       -51.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -36.D0,      0.D0,     0.D0,       20.D0,    0.D0,    0.D0,
     :        44.D0,      0.D0,     0.D0,      -19.D0,    0.D0,    0.D0,
     :       -60.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :        35.D0,      0.D0,     0.D0,      -18.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :        36.D0,      0.D0,     0.D0,      -15.D0,    0.D0,    0.D0,
     :       -36.D0,      0.D0,     0.D0,       20.D0,    0.D0,    0.D0,
     :       -35.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :       -37.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :        35.D0,      0.D0,     0.D0,      -14.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :        65.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :        37.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :       -30.D0,      0.D0,     0.D0,       15.D0,    0.D0,    0.D0,
     :       -32.D0,      0.D0,     0.D0,       16.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=301,320 ) /
     :       -31.D0,      0.D0,     0.D0,       13.D0,    0.D0,    0.D0,
     :        37.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :        31.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :        49.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :       -43.D0,      0.D0,     0.D0,       18.D0,    0.D0,    0.D0,
     :       -32.D0,      0.D0,     0.D0,       14.D0,    0.D0,    0.D0,
     :        30.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -34.D0,      0.D0,     0.D0,       15.D0,    0.D0,    0.D0,
     :       -36.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -38.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -31.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -34.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -35.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        30.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :         0.D0,      0.D0, -1988.D0,        0.D0,    0.D0,-1679.D0,
     :         0.D0,      0.D0,   -63.D0,        0.D0,    0.D0,  -27.D0,
     :         0.D0,      0.D0,   364.D0,        0.D0,    0.D0,  176.D0,
     :         0.D0,      0.D0, -1044.D0,        0.D0,    0.D0, -891.D0,
     :         0.D0,      0.D0,   330.D0,        0.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=321,323 ) /
     :         0.D0,      0.D0,    30.D0,        0.D0,    0.D0,   14.D0,
     :         0.D0,      0.D0,  -162.D0,        0.D0,    0.D0, -138.D0,
     :         0.D0,      0.D0,    75.D0,        0.D0,    0.D0,    0.D0/
 
*
*  Planetary argument multipliers:
*              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre
*
      DATA ( ( NAPL(I,J), I=1,14), J=  1, 20 ) /
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2,
     :         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1,
     :         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0/
      DATA ( ( NAPL(I,J), I=1,14), J= 21, 40 ) /
     :         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
     :         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0/
      DATA ( ( NAPL(I,J), I=1,14), J= 41, 60 ) /
     :         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1/
      DATA ( ( NAPL(I,J), I=1,14), J= 61, 80 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0/
      DATA ( ( NAPL(I,J), I=1,14), J= 81,100 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2/
      DATA ( ( NAPL(I,J), I=1,14), J=101,120 ) /
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1/
      DATA ( ( NAPL(I,J), I=1,14), J=121,140 ) /
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2/
      DATA ( ( NAPL(I,J), I=1,14), J=141,160 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2/
      DATA ( ( NAPL(I,J), I=1,14), J=161,165 ) /
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0/
 
*
*  Planetary nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, cos), obliquity (sin, cos)
*
      DATA ( ( CPL(I,J), I=1,4 ), J=  1, 20 ) /
     :    1440.D0,       0.D0,       0.D0,       0.D0,
     :      56.D0,    -117.D0,     -42.D0,     -40.D0,
     :     125.D0,     -43.D0,       0.D0,     -54.D0,
     :    -114.D0,       0.D0,       0.D0,      61.D0,
     :    -219.D0,      89.D0,       0.D0,       0.D0,
     :    -462.D0,    1604.D0,       0.D0,       0.D0,
     :      99.D0,       0.D0,       0.D0,     -53.D0,
     :      14.D0,    -218.D0,     117.D0,       8.D0,
     :      31.D0,    -481.D0,    -257.D0,     -17.D0,
     :    -491.D0,     128.D0,       0.D0,       0.D0,
     :   -3084.D0,    5123.D0,    2735.D0,    1647.D0,
     :   -1444.D0,    2409.D0,   -1286.D0,    -771.D0,
     :     103.D0,     -60.D0,       0.D0,       0.D0,
     :     -26.D0,     -29.D0,     -16.D0,      14.D0,
     :     284.D0,       0.D0,       0.D0,    -151.D0,
     :     226.D0,     101.D0,       0.D0,       0.D0,
     :     -41.D0,     175.D0,      76.D0,      17.D0,
     :     425.D0,     212.D0,    -133.D0,     269.D0,
     :    1200.D0,     598.D0,     319.D0,    -641.D0,
     :     235.D0,     334.D0,       0.D0,       0.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 21, 40 ) /
     :     266.D0,     -78.D0,       0.D0,       0.D0,
     :    -460.D0,    -435.D0,    -232.D0,     246.D0,
     :       0.D0,     131.D0,       0.D0,       0.D0,
     :     -42.D0,      20.D0,       0.D0,       0.D0,
     :     -10.D0,     233.D0,       0.D0,       0.D0,
     :      78.D0,     -18.D0,       0.D0,       0.D0,
     :      45.D0,     -22.D0,       0.D0,       0.D0,
     :      89.D0,     -16.D0,      -9.D0,     -48.D0,
     :    -349.D0,     -62.D0,       0.D0,       0.D0,
     :     -53.D0,       0.D0,       0.D0,       0.D0,
     :     -21.D0,     -78.D0,       0.D0,       0.D0,
     :      20.D0,     -70.D0,     -37.D0,     -11.D0,
     :      32.D0,      15.D0,      -8.D0,      17.D0,
     :     174.D0,      84.D0,      45.D0,     -93.D0,
     :      11.D0,      56.D0,       0.D0,       0.D0,
     :     -66.D0,     -12.D0,      -6.D0,      35.D0,
     :      47.D0,       8.D0,       4.D0,     -25.D0,
     :      46.D0,      66.D0,      35.D0,     -25.D0,
     :     -68.D0,     -34.D0,     -18.D0,      36.D0,
     :      76.D0,      17.D0,       9.D0,     -41.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 41, 60 ) /
     :      84.D0,     298.D0,     159.D0,     -45.D0,
     :     -82.D0,     292.D0,     156.D0,      44.D0,
     :     -73.D0,      17.D0,       9.D0,      39.D0,
     :    -439.D0,       0.D0,       0.D0,       0.D0,
     :      57.D0,     -28.D0,     -15.D0,     -30.D0,
     :     -40.D0,      57.D0,      30.D0,      21.D0,
     :     273.D0,      80.D0,      43.D0,    -146.D0,
     :    -449.D0,     430.D0,       0.D0,       0.D0,
     :      -8.D0,     -47.D0,     -25.D0,       4.D0,
     :       6.D0,      47.D0,      25.D0,      -3.D0,
     :     -48.D0,    -110.D0,     -59.D0,      26.D0,
     :      51.D0,     114.D0,      61.D0,     -27.D0,
     :    -133.D0,       0.D0,       0.D0,      57.D0,
     :     -18.D0,    -436.D0,    -233.D0,       9.D0,
     :      35.D0,      -7.D0,       0.D0,       0.D0,
     :     -53.D0,      -9.D0,      -5.D0,      28.D0,
     :     -50.D0,     194.D0,     103.D0,      27.D0,
     :     -13.D0,      52.D0,      28.D0,       7.D0,
     :     -91.D0,     248.D0,       0.D0,       0.D0,
     :       6.D0,      49.D0,      26.D0,      -3.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 61, 80 ) /
     :      -6.D0,     -47.D0,     -25.D0,       3.D0,
     :      52.D0,      23.D0,      10.D0,     -23.D0,
     :    -138.D0,       0.D0,       0.D0,       0.D0,
     :      54.D0,       0.D0,       0.D0,     -29.D0,
     :     -37.D0,      35.D0,      19.D0,      20.D0,
     :    -145.D0,      47.D0,       0.D0,       0.D0,
     :     -10.D0,      40.D0,      21.D0,       5.D0,
     :      11.D0,     -49.D0,     -26.D0,      -7.D0,
     :   -2150.D0,       0.D0,       0.D0,     932.D0,
     :      85.D0,       0.D0,       0.D0,     -37.D0,
     :     -86.D0,     153.D0,       0.D0,       0.D0,
     :     -51.D0,       0.D0,       0.D0,      22.D0,
     :     -11.D0,    -268.D0,    -116.D0,       5.D0,
     :      31.D0,       6.D0,       3.D0,     -17.D0,
     :     140.D0,      27.D0,      14.D0,     -75.D0,
     :      57.D0,      11.D0,       6.D0,     -30.D0,
     :     -14.D0,     -39.D0,       0.D0,       0.D0,
     :     -25.D0,      22.D0,       0.D0,       0.D0,
     :      42.D0,     223.D0,     119.D0,     -22.D0,
     :     -27.D0,    -143.D0,     -77.D0,      14.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 81,100 ) /
     :       9.D0,      49.D0,      26.D0,      -5.D0,
     :   -1166.D0,       0.D0,       0.D0,     505.D0,
     :     117.D0,       0.D0,       0.D0,     -63.D0,
     :       0.D0,      31.D0,       0.D0,       0.D0,
     :       0.D0,     -32.D0,     -17.D0,       0.D0,
     :      50.D0,       0.D0,       0.D0,     -27.D0,
     :      30.D0,      -3.D0,      -2.D0,     -16.D0,
     :       8.D0,     614.D0,       0.D0,       0.D0,
     :    -127.D0,      21.D0,       9.D0,      55.D0,
     :     -20.D0,      34.D0,       0.D0,       0.D0,
     :      22.D0,     -87.D0,       0.D0,       0.D0,
     :     -68.D0,      39.D0,       0.D0,       0.D0,
     :       3.D0,      66.D0,      29.D0,      -1.D0,
     :     490.D0,       0.D0,       0.D0,    -213.D0,
     :     -22.D0,      93.D0,      49.D0,      12.D0,
     :     -46.D0,      14.D0,       0.D0,       0.D0,
     :      25.D0,     106.D0,      57.D0,     -13.D0,
     :    1485.D0,       0.D0,       0.D0,       0.D0,
     :      -7.D0,     -32.D0,     -17.D0,       4.D0,
     :      30.D0,      -6.D0,      -2.D0,     -13.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=101,120 ) /
     :     118.D0,       0.D0,       0.D0,     -52.D0,
     :     -28.D0,      36.D0,       0.D0,       0.D0,
     :      14.D0,     -59.D0,     -31.D0,      -8.D0,
     :    -458.D0,       0.D0,       0.D0,     198.D0,
     :       0.D0,     -45.D0,     -20.D0,       0.D0,
     :    -166.D0,     269.D0,       0.D0,       0.D0,
     :     -78.D0,      45.D0,       0.D0,       0.D0,
     :      -5.D0,     328.D0,       0.D0,       0.D0,
     :   -1223.D0,     -26.D0,       0.D0,       0.D0,
     :    -368.D0,       0.D0,       0.D0,       0.D0,
     :     -75.D0,       0.D0,       0.D0,       0.D0,
     :     -13.D0,     -30.D0,       0.D0,       0.D0,
     :     -74.D0,       0.D0,       0.D0,      32.D0,
     :    -262.D0,       0.D0,       0.D0,     114.D0,
     :     202.D0,       0.D0,       0.D0,     -87.D0,
     :      -8.D0,      35.D0,      19.D0,       5.D0,
     :     -35.D0,     -48.D0,     -21.D0,      15.D0,
     :      12.D0,      55.D0,      29.D0,      -6.D0,
     :    -598.D0,       0.D0,       0.D0,       0.D0,
     :       8.D0,     -31.D0,     -16.D0,      -4.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=121,140 ) /
     :     113.D0,       0.D0,       0.D0,     -49.D0,
     :      83.D0,      15.D0,       0.D0,       0.D0,
     :       0.D0,    -114.D0,     -49.D0,       0.D0,
     :     117.D0,       0.D0,       0.D0,     -51.D0,
     :     393.D0,       3.D0,       0.D0,       0.D0,
     :      18.D0,     -29.D0,     -13.D0,      -8.D0,
     :       8.D0,      34.D0,      18.D0,      -4.D0,
     :      89.D0,       0.D0,       0.D0,       0.D0,
     :      54.D0,     -15.D0,      -7.D0,     -24.D0,
     :       0.D0,      35.D0,       0.D0,       0.D0,
     :    -154.D0,     -30.D0,     -13.D0,      67.D0,
     :      80.D0,     -71.D0,     -31.D0,     -35.D0,
     :      61.D0,     -96.D0,     -42.D0,     -27.D0,
     :     123.D0,    -415.D0,    -180.D0,     -53.D0,
     :       0.D0,       0.D0,       0.D0,     -35.D0,
     :       7.D0,     -32.D0,     -17.D0,      -4.D0,
     :     -89.D0,       0.D0,       0.D0,      38.D0,
     :       0.D0,     -86.D0,     -19.D0,      -6.D0,
     :    -123.D0,    -416.D0,    -180.D0,      53.D0,
     :     -62.D0,     -97.D0,     -42.D0,      27.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=141,160 ) /
     :     -85.D0,     -70.D0,     -31.D0,      37.D0,
     :     163.D0,     -12.D0,      -5.D0,     -72.D0,
     :     -63.D0,     -16.D0,      -7.D0,      28.D0,
     :     -21.D0,     -32.D0,     -14.D0,       9.D0,
     :       5.D0,    -173.D0,     -75.D0,      -2.D0,
     :      74.D0,       0.D0,       0.D0,     -32.D0,
     :      83.D0,       0.D0,       0.D0,       0.D0,
     :    -339.D0,       0.D0,       0.D0,     147.D0,
     :      67.D0,     -91.D0,     -39.D0,     -29.D0,
     :      30.D0,     -18.D0,      -8.D0,     -13.D0,
     :       0.D0,    -114.D0,     -50.D0,       0.D0,
     :     517.D0,      16.D0,       7.D0,    -224.D0,
     :     143.D0,      -3.D0,      -1.D0,     -62.D0,
     :      50.D0,       0.D0,       0.D0,     -22.D0,
     :      59.D0,       0.D0,       0.D0,       0.D0,
     :     370.D0,      -8.D0,       0.D0,    -160.D0,
     :      34.D0,       0.D0,       0.D0,     -15.D0,
     :     -37.D0,      -7.D0,      -3.D0,      16.D0,
     :      40.D0,       0.D0,       0.D0,       0.D0,
     :    -184.D0,      -3.D0,      -1.D0,      80.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=161,165 ) /
     :      31.D0,      -6.D0,       0.D0,     -13.D0,
     :      -3.D0,     -32.D0,     -14.D0,       1.D0,
     :     -34.D0,       0.D0,       0.D0,       0.D0,
     :     126.D0,     -63.D0,     -27.D0,     -55.D0,
     :    -126.D0,     -63.D0,     -27.D0,      55.D0/
 
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC
 
*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------
 
*
*  Fundamental (Delaunay) arguments from Simon et al. (1994)
*
      CALL FUNARG ( T,   EL, ELP, F, D, OM )
 
*  Initialize the nutation values.
      DP = 0.D0
      DE = 0.D0
 
*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1
 
*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)
 
*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG
 
 100  CONTINUE
 
*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R
 
*  ------------------
*  PLANETARY NUTATION
*  ------------------
 
*  Planetary longitudes, Mercury through Neptune, wrt mean dynamical
*  ecliptic and equinox of J2000, with high order terms omitted
*  (Simon et al. 1994, 5.8.1-5.8.8).
      ALME = MOD ( 4.402608842461D0 + 2608.790314157421D0 * T, D2PI )
      ALVE = MOD ( 3.176146696956D0 + 1021.328554621099D0 * T, D2PI )
      ALEA = MOD ( 1.753470459496D0 +  628.307584999142D0 * T, D2PI )
      ALMA = MOD ( 6.203476112911D0 +  334.061242669982D0 * T, D2PI )
      ALJU = MOD ( 0.599547105074D0 +   52.969096264064D0 * T, D2PI )
      ALSA = MOD ( 0.874016284019D0 +   21.329910496032D0 * T, D2PI )
      ALUR = MOD ( 5.481293871537D0 +    7.478159856729D0 * T, D2PI )
      ALNE = MOD ( 5.311886286677D0 +    3.813303563778D0 * T, D2PI )
 
*  General precession in longitude (Simon et al. 1994), equivalent
*  to 5028.8200 arcsec/cy at J2000.
      APA = ( 0.024380407358D0 + 0.000005391235D0 * T ) * T
 
*  Initialize the nutation values.
      DP = 0.D0
      DE = 0.D0
 
*  Summation of planetary nutation series (in reverse order).
      DO 200 I = NPL, 1, -1
 
*     Argument and functions.
         ARG = MOD ( DBLE ( NAPL( 1,I) ) * EL   +
     :               DBLE ( NAPL( 2,I) ) * ELP  +
     :               DBLE ( NAPL( 3,I) ) * F    +
     :               DBLE ( NAPL( 4,I) ) * D    +
     :               DBLE ( NAPL( 5,I) ) * OM   +
     :               DBLE ( NAPL( 6,I) ) * ALME +
     :               DBLE ( NAPL( 7,I) ) * ALVE +
     :               DBLE ( NAPL( 8,I) ) * ALEA +
     :               DBLE ( NAPL( 9,I) ) * ALMA +
     :               DBLE ( NAPL(10,I) ) * ALJU +
     :               DBLE ( NAPL(11,I) ) * ALSA +
     :               DBLE ( NAPL(12,I) ) * ALUR +
     :               DBLE ( NAPL(13,I) ) * ALNE +
     :               DBLE ( NAPL(14,I) ) * APA, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)
 
*     Term.
         DP = DP + CPL(1,I) * SARG + CPL(2,I) * CARG
         DE = DE + CPL(3,I) * SARG + CPL(4,I) * CARG
 
 200  CONTINUE
 
*  Convert from 0.1 microarcsec units to radians.
      DPSIPL = DP * U2R
      DEPSPL = DE * U2R
 
*  -----
*  TOTAL
*  -----
 
*  Add planetary and luni-solar components.
      DPSI = DPSIPL + DPSILS
      DEPS = DEPSPL + DEPSLS
 
      RETURN
 
      END
 
 
 
      DOUBLE PRECISION FUNCTION EECT2000 ( DATE1, DATE2 )
*+
*  - - - - - - - - -
*   E E C T 2 0 0 0
*  - - - - - - - - -
*
*  Equation of the equinoxes complementary terms, consistent with
*  IAU 2000 resolutions.
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Capitaine, N., Wallace, P.T., & McCarthy, D.D. (2003). Astron. &
*    Astrophys. 406, pp. 1135-1149, Table 3.
*  IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
*    (Table 5.2e presented in the printed publication is a truncated
*    series. The full series, which is used in NOVAS, is available on
*    the IERS Conventions Center website in file tab5.2e.txt.)
*    ftp://tai.bipm.org/iers/conv2010/chapter5/
*
*  Given:
*     DATE1,DATE2   d    TT date (JD = DATE1+DATE2)
*
*  Returned:
*     EECT00        d    complementary terms (radians)
*
*  This revision:  2002 November 13
*                  References updated 2010 November 26
*
*-----------------------------------------------------------------------
 
      IMPLICIT NONE
 
      DOUBLE PRECISION DATE1, DATE2
 
*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )
 
*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )
 
*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )
 
*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )
 
*  Time since J2000, in Julian centuries
      DOUBLE PRECISION T
 
*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION A, S0, S1
      DOUBLE PRECISION ANMP
 
*  Fundamental arguments
      DOUBLE PRECISION FA(14)
 
*  -----------------------------------------
*  The series for the EE complementary terms
*  -----------------------------------------
 
*  Number of terms in the series
      INTEGER NE0, NE1
      PARAMETER ( NE0=  33, NE1=  1 )
 
*  Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
      INTEGER KE0 ( 14, NE0 ),
     :        KE1 ( 14, NE1 )
 
*  Sine and cosine coefficients
      DOUBLE PRECISION SE0 ( 2, NE0 ),
     :                 SE1 ( 2, NE1 )
 
*  Argument coefficients for t^0
      DATA ( ( KE0(I,J), I=1,14), J =    1,   10 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   11,   20 ) /
     :  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   21,   30 ) /
     :  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1,
     :  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   31,  NE0 ) /
     :  0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
 
*  Argument coefficients for t^1
      DATA ( ( KE1(I,J), I=1,14), J =    1,  NE1 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
 
*  Sine and cosine coefficients for t^0
      DATA ( ( SE0(I,J), I=1,2), J =    1,   10 ) /
     :            +2640.96D-6,          -0.39D-6,
     :              +63.52D-6,          -0.02D-6,
     :              +11.75D-6,          +0.01D-6,
     :              +11.21D-6,          +0.01D-6,
     :               -4.55D-6,          +0.00D-6,
     :               +2.02D-6,          +0.00D-6,
     :               +1.98D-6,          +0.00D-6,
     :               -1.72D-6,          +0.00D-6,
     :               -1.41D-6,          -0.01D-6,
     :               -1.26D-6,          -0.01D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   11,   20 ) /
     :               -0.63D-6,          +0.00D-6,
     :               -0.63D-6,          +0.00D-6,
     :               +0.46D-6,          +0.00D-6,
     :               +0.45D-6,          +0.00D-6,
     :               +0.36D-6,          +0.00D-6,
     :               -0.24D-6,          -0.12D-6,
     :               +0.32D-6,          +0.00D-6,
     :               +0.28D-6,          +0.00D-6,
     :               +0.27D-6,          +0.00D-6,
     :               +0.26D-6,          +0.00D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   21,   30 ) /
     :               -0.21D-6,          +0.00D-6,
     :               +0.19D-6,          +0.00D-6,
     :               +0.18D-6,          +0.00D-6,
     :               -0.10D-6,          +0.05D-6,
     :               +0.15D-6,          +0.00D-6,
     :               -0.14D-6,          +0.00D-6,
     :               +0.14D-6,          +0.00D-6,
     :               -0.14D-6,          +0.00D-6,
     :               +0.14D-6,          +0.00D-6,
     :               +0.13D-6,          +0.00D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   31,  NE0 ) /
     :               -0.11D-6,          +0.00D-6,
     :               +0.11D-6,          +0.00D-6,
     :               +0.11D-6,          +0.00D-6 /
 
*  Sine and cosine coefficients for t^1
      DATA ( ( SE1(I,J), I=1,2), J =    1,  NE1 ) /
     :               -0.87D-6,          +0.00D-6 /
 
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC
 
*  Fundamental Arguments (from IERS Conventions 2000)
 
*  Mean Anomaly of the Moon.
      FA(1) = ANMP ( ( 485868.249036D0 +
     :               ( 715923.2178D0 +
     :               (     31.8792D0 +
     :               (      0.051635D0 +
     :               (     -0.00024470D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1325D0*T, 1D0 ) * D2PI )
 
*  Mean Anomaly of the Sun.
      FA(2) = ANMP ( ( 1287104.793048D0 +
     :               ( 1292581.0481D0 +
     :               (      -0.5532D0 +
     :               (      +0.000136D0 +
     :               (      -0.00001149D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 99D0*T, 1D0 ) * D2PI )
 
*  Mean Longitude of the Moon minus Mean Longitude of the Ascending
*  Node of the Moon.
      FA(3) = ANMP ( (  335779.526232D0 +
     :               (  295262.8478D0 +
     :               (     -12.7512D0 +
     :               (      -0.001037D0 +
     :               (       0.00000417D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1342D0*T, 1D0 ) * D2PI )
 
*  Mean Elongation of the Moon from the Sun.
      FA(4) = ANMP ( ( 1072260.703692D0 +
     :               ( 1105601.2090D0 +
     :               (      -6.3706D0 +
     :               (       0.006593D0 +
     :               (      -0.00003169D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1236D0*T, 1D0 ) * D2PI )
 
*  Mean Longitude of the Ascending Node of the Moon.
      FA(5) = ANMP ( (  450160.398036D0 +
     :               ( -482890.5431D0 +
     :               (       7.4722D0 +
     :               (       0.007702D0 +
     :               (      -0.00005939D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( -5D0*T, 1D0 ) * D2PI )
 
      FA( 6) = ANMP ( 4.402608842D0 + 2608.7903141574D0 * T )
      FA( 7) = ANMP ( 3.176146697D0 + 1021.3285546211D0 * T )
      FA( 8) = ANMP ( 1.753470314D0 +  628.3075849991D0 * T )
      FA( 9) = ANMP ( 6.203480913D0 +  334.0612426700D0 * T )
      FA(10) = ANMP ( 0.599546497D0 +   52.9690962641D0 * T )
      FA(11) = ANMP ( 0.874016757D0 +   21.3299104960D0 * T )
      FA(12) = ANMP ( 5.481293872D0 +    7.4781598567D0 * T )
      FA(13) = ANMP ( 5.311886287D0 +    3.8133035638D0 * T )
      FA(14) =      ( 0.024381750D0 +    0.00000538691D0 * T ) * T
 
*  Evaluate the EE complementary terms.
      S0 = 0D0
      S1 = 0D0
 
      DO I = NE0,1,-1
         A = 0D0
         DO J=1,14
            A = A + DBLE(KE0(J,I))*FA(J)
         END DO
         S0 = S0 + ( SE0(1,I)*SIN(A) + SE0(2,I)*COS(A) )
      END DO
      DO I = NE1,1,-1
         A = 0D0
         DO J=1,14
            A = A + DBLE(KE1(J,I))*FA(J)
         END DO
         S1 = S1 + ( SE1(1,I)*SIN(A) + SE1(2,I)*COS(A) )
      END DO
      EECT2000 = ( S0 + S1 * T ) * DAS2R
 
*  Finished.
 
      END
 
 
 
      DOUBLE PRECISION FUNCTION ANMP ( A )
 
*  Normalize angle into the range -pi <= A < +pi.
 
      IMPLICIT NONE
 
      DOUBLE PRECISION A
 
      DOUBLE PRECISION DPI, D2PI
      PARAMETER ( DPI = 3.141592653589793238462643D0,
     :            D2PI = 6.283185307179586476925287D0 )
 
      DOUBLE PRECISION W
 
      W = MOD(A,D2PI)
      IF ( ABS(W) .GE. DPI ) W = W - SIGN(D2PI,A)
      ANMP = W
 
      END
 
c*********************************************
c****
c**** #include "NOVAS_F3.1_solsys3.f"
c****
c*********************************************
 
*** SOLSYS VERSION 3 PACKAGE: SOLSYS, SUN, IDSS ***
 
      SUBROUTINE SOLSYS (TJD,M,K,POS,VEL,IERR)
*
*     SUBROUTINE SOLSYS VERSION 3.
*     THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF THE
*     EARTH AT EPOCH TJD BY EVALUATING A CLOSED-FORM THEORY WITHOUT
*     REFERENCE TO AN EXTERNAL FILE.  THIS ROUTINE CAN ALSO PROVIDE
*     THE POSITION AND VELOCITY OF THE SUN.
*
*          TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
*          M    = BODY IDENTIFICATION NUMBER (IN)
*                 SET M=0 OR M=1 FOR THE SUN
*                 SET M=2 OR M=3 FOR THE EARTH
*          K    = ORIGIN SELECTION CODE (IN)
*                 SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
*                 SET K=1 FOR ORIGIN AT CENTER OF SUN
*          POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
*                 OF J2000.0, COMPONENTS IN AU (OUT)
*          VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
*                 OF J2000.0, COMPONENTS IN AU/DAY (OUT)
*          IERR = ERROR INDICATOR (OUT)
*                 IERR=0 MEANS EVERYTHING OK
*                 IERR=1 MEANS TJD BEFORE FIRST ALLOWED DATE
*                 IERR=2 MEANS TJD AFTER LAST ALLOWED DATE
*
*
       DOUBLE PRECISION TJD,POS,VEL,PI,TWOPI,T0,OBL,EL,C,P,TLAST,
     .     PM,PA,PE,PJ,PO,PW,PL,PN,
     .     TMASS,SE,CE,SI,CI,SN,CN,SW,CW,P1,P2,P3,Q1,Q2,Q3,ROOTE,A,B,
     .     QJD,E,MLON,MA,U,SINU,COSU,ANR,PPLAN,VPLAN,F,PBARY,VBARY,
     .     DFLOAT,DABS,DMOD,DSIN,DCOS,DSQRT
      DIMENSION POS(3), VEL(3), EL(21), C(13), P(3,3),
     .     PM(4), PA(4), PE(4), PJ(4), PO(4), PW(4), PL(4), PN(4),
     .     A(3,4), B(3,4), PPLAN(3), VPLAN(3), PBARY(3), VBARY(3)
      SAVE
 
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( TWOPI  = 2.D0 * PI             )
      PARAMETER ( T0     = 2451545.0D0           )
      PARAMETER ( OBL    = 23.43927944D0         )
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0
*     OBL = OBLIQUITY OF ECLIPTIC AT EPOCH J2000.0
 
      DATA EL, C, P / 43*0.D0 /,   TLAST / 0.D0 /
 
*     ARRAYS BELOW CONTAIN MASSES AND ORBITAL ELEMENTS OF THE FOUR
*     LARGEST PLANETS (SEE EXPLANATORY SUPPLEMENT (1992), P. 316)
*     WITH ANGLES IN RADIANS
*     THIS DATA USED FOR BARYCENTER COMPUTATIONS ONLY
*                 JUPITER        SATURN        URANUS       NEPTUNE
      DATA PM /  1047.349D 0,  3497.898D 0,   22903.0D 0,   19412.2D 0 /
      DATA PA /  5.203363D 0,  9.537070D 0, 19.191264D 0, 30.068963D 0 /
      DATA PE /  0.048393D 0,  0.054151D 0,  0.047168D 0,  0.008586D 0 /
      DATA PJ /  0.022782D 0,  0.043362D 0,  0.013437D 0,  0.030878D 0 /
      DATA PO /  1.755036D 0,  1.984702D 0,  1.295556D 0,  2.298977D 0 /
      DATA PW /  0.257503D 0,  1.613242D 0,  2.983889D 0,  0.784898D 0 /
      DATA PL /  0.600470D 0,  0.871693D 0,  5.466933D 0,  5.321160D 0 /
      DATA PN /  1.450138D-3,  5.841727D-4,  2.047497D-4,  1.043891D-4 /
 
      IF ( TLAST .LT. 1.D0 ) THEN
*         FIRST TIME COMPUTATIONS
*         MASS OF SUN PLUS FOUR INNER PLANETS
          TMASS = 1.D0 + 5.977D-6
          SE = DSIN ( OBL * PI / 180.D0 )
          CE = DCOS ( OBL * PI / 180.D0 )
          DO 15 I = 1, 4
              TMASS = TMASS + 1.D0 / PM(I)
*             COMPUTE SINE AND COSINE OF ORBITAL ANGLES
			  SI = DSIN ( PJ(I) )
			  CI = DCOS ( PJ(I) )
			  SN = DSIN ( PO(I) )
			  CN = DCOS ( PO(I) )
			  SW = DSIN ( PW(I) - PO(I) )
			  CW = DCOS ( PW(I) - PO(I) )
*             COMPUTE P AND Q VECTORS (SEE BROUWER & CLEMENCE (1961),
*             METHODS OF CELESTIAL MECHANICS, PP. 35-36.)
			  P1 =    CW * CN - SW * SN * CI
			  P2 = (  CW * SN + SW * CN * CI ) * CE - SW * SI * SE
			  P3 = (  CW * SN + SW * CN * CI ) * SE + SW * SI * CE
			  Q1 =   -SW * CN - CW * SN * CI
			  Q2 = ( -SW * SN + CW * CN * CI ) * CE - CW * SI * SE
			  Q3 = ( -SW * SN + CW * CN * CI ) * SE + CW * SI * CE
              ROOTE = DSQRT ( 1.D0 - PE(I)**2 )
              A(1,I) = PA(I) * P1
              A(2,I) = PA(I) * P2
              A(3,I) = PA(I) * P3
              B(1,I) = PA(I) * ROOTE * Q1
              B(2,I) = PA(I) * ROOTE * Q2
              B(3,I) = PA(I) * ROOTE * Q3
  15      CONTINUE
          TLAST = 1.D0
      END IF
 
      IERR = 0
*     VALID DATES ARE WITHIN 3 CENTURIES OF J2000, ALTHOUGH RESULTS
*     DETERIORATE GRADUALLY
      IF ( TJD .LT. 2340000.5D0 ) IERR = 1
      IF ( TJD .GT. 2560000.5D0 ) IERR = 2
      IF ( IERR .NE. 0 ) GO TO 110
      IF ( M .GE. 2 ) GO TO 30
 
*     FORM HELIOCENTRIC COORDINATES OF SUN
  20  DO 25 J=1,3
          POS(J) = 0.D0
          VEL(J) = 0.D0
  25  CONTINUE
      IF ( K .GE. 1 ) GO TO 110
      GO TO 90
 
*     FORM HELIOCENTRIC COORDINATES OF EARTH
*     VELOCITIES ARE OBTAINED FROM CRUDE NUMERICAL DIFFERENTIATION
  30  DO 35 I = 1, 3
          QJD = TJD + DFLOAT(I-2) * 0.1D0
C         SUBROUTINE SUN COMPUTES EARTH-SUN VECTOR
          CALL SUN ( QJD, EL, C )
          CALL PRECES ( QJD, C(11), T0, POS )
          P(I,1) = -POS(1)
          P(I,2) = -POS(2)
          P(I,3) = -POS(3)
  35  CONTINUE
      DO 40 J=1,3
          POS(J) =   P(2,J)
          VEL(J) = ( P(3,J) - P(1,J) ) / 0.2D0
  40  CONTINUE
      IF ( K .GE. 1 ) GO TO 110
 
*     IF K=0, MOVE ORIGIN TO SOLAR SYSTEM BARYCENTER
*     SOLAR SYSTEM BARYCENTER COORDINATES COMPUTED FROM KEPLERIAN
*     APPROXIMATIONS OF THE COORDINATES OF THE FOUR LARGEST PLANETS
  90  IF ( DABS ( TJD - TLAST ) .LT. 1.D-6 ) GO TO 99
      DO 92 J = 1, 3
          PBARY(J) = 0.D0
          VBARY(J) = 0.D0
  92  CONTINUE
*     THE FOLLOWING LOOP CYCLES ONCE FOR EACH OF THE FOUR LARGE PLANETS
      DO 98 I = 1, 4
*         COMPUTE MEAN LONGITUDE, MEAN ANOMALY, AND ECCENTRIC ANOMOLY
          E = PE(I)
          MLON = PL(I) + PN(I) * ( TJD - T0 )
          MA = DMOD ( MLON - PW(I), TWOPI )
          U = MA + E * DSIN ( MA ) + 0.5D0 * E * E * DSIN ( 2.D0 * MA )
          SINU = DSIN ( U )
          COSU = DCOS ( U )
*         COMPUTE VELOCITY FACTOR
          ANR = PN(I) / ( 1.D0 - E * COSU )
*         COMPUTE PLANET'S POSITION AND VELOCITY WRT EQ & EQ J2000
          PPLAN(1) = A(1,I) * ( COSU - E ) + B(1,I) * SINU
          PPLAN(2) = A(2,I) * ( COSU - E ) + B(2,I) * SINU
          PPLAN(3) = A(3,I) * ( COSU - E ) + B(3,I) * SINU
          VPLAN(1) = ANR * ( -A(1,I) * SINU + B(1,I) * COSU )
          VPLAN(2) = ANR * ( -A(2,I) * SINU + B(2,I) * COSU )
          VPLAN(3) = ANR * ( -A(3,I) * SINU + B(3,I) * COSU )
*         COMPUTE MASS FACTOR AND ADD IN TO TOTAL DISPLACEMENT
          F = 1.D0 / ( PM(I) * TMASS )
          PBARY(1) = PBARY(1) + PPLAN(1) * F
          PBARY(2) = PBARY(2) + PPLAN(2) * F
          PBARY(3) = PBARY(3) + PPLAN(3) * F
          VBARY(1) = VBARY(1) + VPLAN(1) * F
          VBARY(2) = VBARY(2) + VPLAN(2) * F
          VBARY(3) = VBARY(3) + VPLAN(3) * F
  98  CONTINUE
      TLAST = TJD
  99  DO 100 J=1,3
          POS(J) = POS(J) - PBARY(J)
          VEL(J) = VEL(J) - VBARY(J)
 100  CONTINUE
 
 110  RETURN
 
      END
 
 
 
      SUBROUTINE SUN (DJ,EL,C)
C
C     FOR USE WITH SUBROUTINE SOLSYS VERSION 3.
C     THIS SUBROUTINE COMPUTES THE COORDINATES OF THE EARTH-SUN
C     POSITION VECTOR WITH RESPECT TO THE ECLIPTIC AND EQUATOR
C     OF DATE.  A MODIFIED FORM OF NEWCOMB'S THEORY ('TABLES OF THE
C     SUN', 1898) IS USED.  ONLY THE LARGEST PERIODIC PERTURBATIONS
C     ARE EVALUATED, AND VAN FLANDERN'S EXPRESSIONS FOR THE FUNDAMENTAL
C     ARGUMENTS ('IMPROVED MEAN ELEMENTS FOR THE EARTH AND MOON', 1981)
C     ARE USED.  THE ABSOLUTE ACCURACY IS NO WORSE THAN 1 ARCSECOND
C     (AVERAGE ERROR ABOUT 0.2 ARCSECOND) OVER 1800-2200.
C     (ADAPTED FROM SUBROUTINE IAUSUN BY P. M. JANICZEK, USNO.)
C
C          DJ   = TDB JULIAN DATE OF DESIRED EPOCH (IN)
C          EL   = ARRAY OF ORBITAL ELEMENTS (SEE BELOW) FOR
C                 EPOCH DJ (OUT)
C          C    = ARRAY OF COORDINATES (SEE BELOW) FOR
C                 EPOCH DJ (OUT)
C
C
      DOUBLE PRECISION DJ,EL,C,T,TP,T20,RO,GV,GM,GJ,GS,DL,DR,DB,DG,
     1 DBLARG,D,TWOPI,STR,RTD,R,TR,
     2 SINO,COSO,SINL,COSL,SINB,COSB,
     3 DSIN,DCOS,DMOD
C
      DIMENSION EL(21)
C
C     EL( 1)= SEMI-MAJOR AXIS, AU
C     EL( 2)= ORBITAL ECCENTRICITY
C     EL( 5)= LONGITUDE OF PERIGEE, RADIANS
C     EL( 9)= UNPERTURBED MEAN LONGITUDE, RADIANS
C     EL(10)= MEAN ANOMALY, AFFECTED BY LONG-PD PERTURBATIONS, RADIANS
C     EL(11)= UNPERTURBED RADIUS, AU
C     EL(12)= EQUATION OF THE CENTER, RADIANS
C     EL(13)= MEAN OBLIQUITY OF ECLIPTIC, RADIANS
C     EL(14)= MEAN LONGITUDE OF MOON, RADIANS
C     EL(15)= MEAN ANOMALY OF MOON, RADIANS
C     EL(16)= LUNAR MEAN ARGUMENT OF LATITUDE, RADIANS
C     EL(17)= MEAN LONGITUDE OF LUNAR ASCENDING NODE, RADIANS
C     EL(21)= MEAN LONGITUDE OF MOON'S PERIGEE, RADIANS
C             (REMAINING ELEMENTS OF ARRAY EL NOT USED)
C
      DIMENSION C(13)
C
C     C( 1) = PERTURBED RADIUS VECTOR, AU
C     C( 2) = SAME AS C(4), DEGREES
C     C( 3) = SAME AS C(5), DEGREES
C     C( 4) = ECLIPTIC LONGITUDE WRT MEAN ECL & EQUX OF DATE, RADIANS
C     C( 5) = ECLIPTIC LATITUDE  WRT MEAN ECL        OF DATE, RADIANS
C     C(11) = EQUATORIAL X WRT MEAN EQU & EQUX OF DATE, AU
C     C(12) = EQUATORIAL Y WRT MEAN EQU & EQUX OF DATE, AU
C     C(13) = EQUATORIAL Z WRT MEAN EQU & EQUX OF DATE, AU
C             (REMAINING ELEMENTS OF ARRAY C NOT USED)
C
C
C***********************************************************************
C
C     PART I    TABLES OF THE PERTURBATIONS
C
      DIMENSION X(8,46), X1(80), X2(80), X3(80), X4(80), X5(48)
      EQUIVALENCE (X(1, 1),X1(1))
      EQUIVALENCE (X(1,11),X2(1))
      EQUIVALENCE (X(1,21),X3(1))
      EQUIVALENCE (X(1,31),X4(1))
      EQUIVALENCE (X(1,41),X5(1))
C
C     PERTURBATIONS BY VENUS
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA X1 /  - 1.,  0., +  33.,-  67., -  85.,-  39., +  24.,-  17.,
     2           - 1.,+ 1., +2353.,-4228., -2062.,-1146., -   4.,+   3.,
     3           - 1.,+ 2., -  65.,-  34., +  68.,-  14., +   6.,-  92.,
     4           - 2.,+ 1., -  99.,+  60., +  84.,+ 136., +  23.,-   3.,
     5           - 2.,+ 2., -4702.,+2903., +3593.,+5822., +  10.,-   6.,
     6           - 2.,+ 3., +1795.,-1737., - 596.,- 632., +  37.,-  56.,
     7           - 3.,+ 3., - 666.,+  27., +  44.,+1044., +   8.,+   1.,
     8           - 3.,+ 4., +1508.,- 397., - 381.,-1448., + 185.,- 100.,
     9           - 3.,+ 5., + 763.,- 684., + 126.,+ 148., +   6.,-   3.,
     *           - 4.,+ 4., - 188.,-  93., - 166.,+ 337.,     0.,    0./
      DATA X2 /  - 4.,+ 5., - 139.,-  38., -  51.,+ 189., -  31.,-   1.,
     2           - 4.,+ 6., + 146.,-  42., -  25.,-  91., +  12.,    0.,
     3           - 5.,+ 5., -  47.,-  69., - 134.,+  93.,     0.,    0.,
     4           - 5.,+ 7., - 119.,-  33., -  37.,+ 136., -  18.,-   6.,
     5           - 5.,+ 8., + 154.,    0.,     0.,-  26.,     0.,    0.,
     6           - 6.,+ 6., -   4.,-  38., -  80.,+   8.,     0.,    0.,
C
C     PERTURBATIONS BY MARS
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
     7           + 1.,- 1., - 216.,- 167., -  92.,+ 119.,     0.,    0.,
     8           + 2.,- 2., +1963.,- 567., - 573.,-1976.,     0.,-   8.,
     9           + 2.,- 1., -1659.,- 617., +  64.,- 137.,     0.,    0.,
     *           + 3.,- 3., +  53.,- 118., - 154.,-  67.,     0.,    0./
      DATA X3 /  + 3.,- 2., + 396.,- 153., -  77.,- 201.,     0.,    0.,
     2           + 4.,- 3., - 131.,+ 483., + 461.,+ 125., +   7.,+   1.,
     3           + 4.,- 2., + 526.,- 256., +  43.,+  96.,     0.,    0.,
     4           + 5.,- 4., +  49.,+  69., +  87.,-  62.,     0.,    0.,
     5           + 5.,- 3., -  38.,+ 200., +  87.,+  17.,     0.,    0.,
     6           + 6.,- 4., - 104.,- 113., - 102.,+  94.,     0.,    0.,
     7           + 6.,- 3., -  11.,+ 100., -  27.,-   4.,     0.,    0.,
     8           + 7.,- 4., -  78.,-  72., -  26.,+  28.,     0.,    0.,
     9           + 9.,- 5., +  60.,-  15., -   4.,-  17.,     0.,    0.,
     *           +15.,- 8., + 200.,-  30., -   1.,-   6.,     0.,    0./
C
C     PERTURBATIONS BY JUPITER
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA X4 /  + 1.,- 2., - 155.,-  52., -  78.,+ 193., +   7.,    0.,
     2           + 1.,- 1., -7208.,+  59., +  56.,+7067., -   1.,+  17.,
     3           + 1.,  0., - 307.,-2582., + 227.,-  89., +  16.,    0.,
     4           + 1.,+ 1., +   8.,-  73., +  79.,+   9., +   1.,+  23.,
     5           + 2.,- 3., +  11.,+  68., + 102.,-  17.,     0.,    0.,
     6           + 2.,- 2., + 136.,+2728., +4021.,- 203.,     0.,    0.,
     7           + 2.,- 1., - 537.,+1518., +1376.,+ 486., +  13.,+ 166.,
     8           + 3.,- 3., - 162.,+  27., +  43.,+ 278.,     0.,    0.,
     9           + 3.,- 2., +  71.,+ 551., + 796.,- 104., +   6.,-   1.,
     *           + 3.,- 1., -  31.,+ 208., + 172.,+  26., +   1.,+  18./
      DATA X5 /  + 4.,- 3., -  43.,+   9., +  13.,+  73.,     0.,    0.,
     2           + 4.,- 2., +  17.,+  78., + 110.,-  24.,     0.,    0.,
C
C     PERTURBATIONS BY SATURN
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
     3           + 1.,- 1., -  77.,+ 412., + 422.,+  79., +   1.,+   6.,
     4           + 1.,  0., -   3.,- 320., +   8.,-   1.,     0.,    0.,
     5           + 2.,- 2., +  38.,- 101., - 152.,-  57.,     0.,    0.,
     6           + 2.,- 1., +  45.,- 103., - 103.,-  44.,     0.,    0./
C
C
C***********************************************************************
C
C     PART II   NECESSARY PRELIMINARIES
C
      DATA TWOPI /6.283185307179586D0/
      DATA STR   /206264806.2470964D0/
      DATA RTD   /57.295779513082321D0/
      DATA R     /1296000.0D0/
      TR = 1000.0D0 / STR
C
C     T  = TIME IN JULIAN CENTURIES FROM 1900 JANUARY 0
      T  = (DJ - 2415020.D0)/36525.D0
C
C     TP = TIME IN JULIAN YEARS     FROM 1850 JANUARY 0
      TP = (DJ - 2396758.D0)/365.25D0
C
C     T20= TIME IN JULIAN CENTURIES FROM J2000.0
      T20= (DJ - 2451545.D0)/36525.D0
C
C
C***********************************************************************
C
C     PART III  COMPUTATION OF ELLIPTIC ELEMENTS AND SECULAR TERMS
C
C     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
      EL( 1) = 1.00000030007166D0
      EL( 2) = 0.016708320D0 + (-0.42229D-04 - 0.126D-06 * T20) * T20
      EL( 5) = 1018578.046D0 + (6190.046D0 +
     1                (1.666D0 + 0.012D0 * T20) * T20) * T20
      EL( 5) = EL( 5) * TR
      EL( 9) = 1009677.850D0 + (100.0D0 * R + 2771.27D0 +
     1                1.089D0 * T20) * T20
      EL( 9) = DMOD (EL( 9) * TR, TWOPI)
      EL(10) = 1287099.804D0 + (99.0D0 * R + 1292581.224D0 +
     1                (-0.577D0 - 0.012D0 * T20) * T20) * T20
      EL(10) = DMOD (EL(10) * TR, TWOPI)
 
C     EXPRESSION FOR OBLIQUITY FROM P03 (IAU 2006) PRECESSION
      EL(13) = 84381.406D0 + (-46.836769D0 +
     1               (-0.0001831D0 + 0.00200340D0 * T20) * T20) * T20
      EL(13) = EL(13) * TR
 
C     KAPLAN CORRECTION TO SUN'S MEAN LONGITUDE TO FIT DE405 OVER
C     INTERVAL 1800-2200, USING P03 (IAU 2006) PRECESSION
      EL(9) = EL(9)
     1      + ( 0.1320D0 - 0.1355D0 * T20 ) * TR
 
C
C***********************************************************************
C
C     PART IV   LUNAR TERMS
C
C     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
      EL(14) = 785939.157D0 + (1336.0D0 * R + 1108372.598D0
     1                + (-5.802D0 + 0.019D0 * T20) * T20) * T20
      EL(14) = DMOD (EL(14) * TR, TWOPI)
      EL(17) = 450160.280D0 + (-5.0D0 * R - 482890.539D0 +
     1                (7.455D0 + 0.008D0 * T20) * T20) * T20
      EL(17) = DMOD (EL(17) * TR, TWOPI)
      EL(21) = 300072.424D0 + (11.0D0 * R + 392449.965D0 +
     1                (-37.112D0 - 0.045D0 * T20) * T20) * T20
      EL(21) = DMOD (EL(21) * TR, TWOPI)
C
C     DERIVED ARGUMENTS
      EL(15) = EL(14) - EL(21)
      EL(16) = EL(14) - EL(17)
      EL(15) = DMOD (EL(15),TWOPI)
      EL(16) = DMOD (EL(16),TWOPI)
C     MEAN ELONGATION
      D      = EL(14) - EL(9)
C
C     COMBINATIONS OF ARGUMENTS AND THE PERTURBATIONS
      D = DMOD (D,TWOPI)
      ARG = D
      DL =    +  6469.*SIN(ARG) +  13.*SIN(3.*ARG)
      DR =    + 13390.*COS(ARG) +  30.*COS(3.*ARG)
C
      DBLARG = D + EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +  177.*SIN(ARG)
      DR = DR +  370.*COS(ARG)
C
      DBLARG = D - EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL -  424.*SIN(ARG)
      DR = DR - 1330.*COS(ARG)
C
      DBLARG = 3.D0*D - EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +   39.*SIN(ARG)
      DR = DR +   80.*COS(ARG)
C
      DBLARG = D + EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL -   64.*SIN(ARG)
      DR = DR -  140.*COS(ARG)
C
      DBLARG = D - EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +  172.*SIN(ARG)
      DR = DR +  360.*COS(ARG)
C
      EL(16) = DMOD (EL(16),TWOPI)
      ARG = EL(16)
      DB =    + 576.*SIN(ARG)
C
C
C***********************************************************************
C
C     PART V    COMPUTATION OF PERIODIC PERTURBATIONS
C
C     THE PERTURBING MEAN ANOMALIES
C
      GV  = 0.19984020D+01 + .1021322923D+02*TP
      GM  = 0.19173489D+01 + .3340556174D+01*TP
      GJ  = 0.25836283D+01 + .5296346478D+00*TP
      GS  = 0.49692316D+01 + .2132432808D+00*TP
      GV  = DMOD (GV,TWOPI)
      GM  = DMOD (GM,TWOPI)
      GJ  = DMOD (GJ,TWOPI)
      GS  = DMOD (GS,TWOPI)
C
C
C     MODIFICATION OF FUNDAMENTAL ARGUMENTS
C
C     APPLICATION OF THE JUPITER-SATURN GREAT INEQUALITY
C     TO JUPITER'S MEAN ANOMALY
C
      GJ = GJ + 0.579904067D-02 * DSIN (5.D0*GS - 2.D0*GJ
     1                 + 1.1719644977D0 - 0.397401726D-03*TP)
      GJ = DMOD (GJ,TWOPI)
C
C     LONG PERIOD PERTURBATIONS OF MEAN ANOMALY
C
      ST = T
C                ARGUMENT IS ( 4 MARS - 7 EARTH + 3 VENUS )
      DG = 266.* SIN (0.555015 + 2.076942*ST)
C                ARGUMENT IS ( 3 JUPITER - 8 MARS + 4 EARTH )
     1    + 6400.* SIN (4.035027 + 0.3525565*ST)
C                ARGUMENT IS ( 13 EARTH - 8 VENUS )
     2    + (1882.-16.*ST) * SIN (0.9990265 + 2.622706*ST)
C
C
C     COMPUTATION OF THE EQUATION OF THE CENTER
C
C     FORM PERTURBED MEAN ANOMALY
      EL(10) = DG/STR + EL(10)
      EL(10) = DMOD (EL(10),TWOPI)
      EL(12) =   DSIN(     EL(10)) * (6910057.D0 -(17240.D0+52.D0*T)*T)
     1         + DSIN(2.D0*EL(10)) * (  72338.D0 -    361.D0*T)
     2         + DSIN(3.D0*EL(10)) * (   1054.D0 -      1.D0*T)
C
C     THE UNPERTURBED RADIUS VECTOR
      RO     =                          30570.D0 -    150.D0*T
     1         - DCOS(     EL(10)) * (7274120.D0 - (18140.D0+50.D0*T)*T)
     2         - DCOS(2.D0*EL(10)) * (  91380.D0 -    460.D0*T)
     3         - DCOS(3.D0*EL(10)) * (   1450.D0 -     10.D0*T)
      EL(11) = 10.D0**(RO*1.D-09)
C
C
C     SELECTED PLANETARY PERTURBATIONS FROM NEWCOMB'S THEORY FOLLOW
C
C     PERTURBATIONS BY VENUS
      DO 20 K=1,16
C     ARGUMENT J * VENUS +   I * EARTH
      DBLARG = X(1,K)*GV + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   20 CONTINUE
C
C     PERTURBATIONS BY MARS
      DO 30 K=17,30
C     ARGUMENT  J * MARS +   I * EARTH
      DBLARG = X(1,K)*GM + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   30 CONTINUE
C
C     PERTURBATIONS BY JUPITER
      DO 40 K=31,42
C     ARGUMENT J*JUPITER +   I * EARTH
      DBLARG = X(1,K)*GJ + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   40 CONTINUE
C
C     PERTURBATIONS BY SATURN
      DO 50 K=43,46
C     ARGUMENT J*SATURN  +   I * EARTH
      DBLARG = X(1,K)*GS + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   50 CONTINUE
C
C
C***********************************************************************
C
C     PART VI   COMPUTATION OF ECLIPTIC AND EQUATORIAL COORDINATES
C
      C(1) = EL(11)*10.D0**(DR*1.D-09)
      C(4) = (DL + DG + EL(12))/STR + EL(9)
      C(4) = DMOD (C(4),TWOPI)
      C(5) = DB/STR
      C(2) = C(4)*RTD
      C(3) = C(5)*RTD
      SINO = DSIN (EL(13))
      COSO = DCOS (EL(13))
      SINL = DSIN (C(4))
      COSL = DCOS (C(4))
      SINB = DSIN (C(5))
      COSB = DCOS (C(5))
      C(11) = C(1) * (COSB * COSL)
      C(12) = C(1) * (COSB * SINL * COSO - SINB * SINO)
      C(13) = C(1) * (COSB * SINL * SINO + SINB * COSO)
C
C
C***********************************************************************
C
C
      RETURN
C
      END
 
 
 
      INTEGER FUNCTION IDSS ( NAME )
*
*     THIS FUNCTION RETURNS THE ID NUMBER OF A SOLAR SYSTEM BODY
*     FOR THE VERSION OF SOLSYS (OR SOLSYS-AUXPOS COMBINATION) IN USE.
*
*         NAME   = NAME OF BODY WHOSE ID NUMBER IS DESIRED, E.G.,
*                  'SUN', 'MOON, 'MERCURY', ETC., EXPRESSED AS ALL
*                  UPPER-CASE LETTERS (IN)
*         IDSS   = ID NUMBER OF BODY, FOR USE IN CALLS TO SOLSYS
*                  (FUNCTION VALUE RETURNED)
*
*     NOTE 1: IN THIS VERSION, ONLY THE FIRST THREE LETTERS OF THE
*     BODY'S NAME ARE USED FOR IDENTIFICATION.  ALTERNATIVE VERSIONS
*     MIGHT USE MORE LETTERS.
*
*     NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=1, SINCE SOLSYS
*     VERSION 3 DOES NOT PROCESS SPLIT JULIAN DATES.
*
*     NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
*     THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.
*
*
      CHARACTER NAME*(*), NAMEIN*3, NAMES*3
      DIMENSION NAMES(35), IDS(35)
 
      DATA NAMES / 'SUN', 'EAR', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---'  /
      DATA IDS   /     0,     3,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0  /
      DATA NUM   / 2 /
 
   3  FORMAT ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', A )
 
      IDSS = -9999
      NAMEIN = NAME
 
*     LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
      DO 20 I = 1, NUM
          IF ( NAMEIN .EQ. NAMES(I) ) THEN
              IDSS = IDS(I)
              GO TO 30
          END IF
  20  CONTINUE
 
*     IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES
      IF ( NAMEIN .EQ. 'JD ' ) THEN
*         IN THIS CASE, SET IDSS=2 IF SOLSYS PROCESSES SPLIT
*         JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE
          IDSS = 1
          GO TO 30
      END IF
 
      WRITE ( *, 3 ) NAME
 
  30  RETURN
 
      END
 
c*********************************************
c****
c**** #include "/user/jayander/FORTRAN/WCS/def_pi.f"
c****
c*********************************************
 
#define _PI_  3.14159265358979323846d0
 
 
c-----------------------------------------------
c
c This routine will take a point (x,y) in the
c tangent plane, and the coordinate of the
c tangent-point (ra0,dec0), and will compute
c the ra coordinate that corresponds to (x,y).
c
      real*8 function xy2r(x,y,r0,d0)
      implicit none
 
      real*8 x, y
      real*8 r0,d0
 
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 tandr, dr
      real*8 xrad, yrad
 
      xrad = x/180.0d0*(_PI_)
      yrad = y/180.0d0*(_PI_)
 
      cosd0 = cos(d0*(_PI_)/180.0d0)
      sind0 = sin(d0*(_PI_)/180.0d0)
 
      tandr = xrad/(cosd0-yrad*sind0)
      dr    = atan(tandr)
      xy2r  = r0 + dr*180.0d0/(_PI_)
 
      return
      end
 
 
c-----------------------------------------------
c
c This routine will take a point (x,y) in the
c tangent plane, and the coordinate of the
c tangent-point (ra0,dec0), and will compute
c the ra coordinate that corresponds to (x,y).
c
      real*8 function xy2d(x,y,r0,d0)
      implicit none
 
      real*8 x, y
      real*8 r0,d0
 
      real*8 cosd0, sind0
      real*8 tandr, dr
      real*8 cosdr
      real*8 tande
      real*8 xrad, yrad
 
      xrad = x/180.0d0*(_PI_)
      yrad = y/180.0d0*(_PI_)
 
      cosd0 = cos(d0*(_PI_)/180.0d0)
      sind0 = sin(d0*(_PI_)/180.0d0)
 
      tandr = xrad/(cosd0-yrad*sind0)
      dr    = atan(tandr)
 
      cosdr = cos(dr)
 
      tande = cosdr*(sind0+yrad*cosd0)/
     .              (cosd0-yrad*sind0)
 
      !print*,'---> tande: ',tande
      xy2d = atan(tande)*180.0d0/(_PI_)
 
      return
      end
 
 
 
 
c---------------------------------
c compact 50-byte storage:
c---------------------------------
c
c  01-07  ra   + error
c  08-14  de   + error
c  15-19  pi   + error
c  20-24  pmra + error
c  25-29  pmde + error
c  30-33  gmag + error
c  34-37  gmag + error
c  38-41  gmag + error
c  42-49  source_id
c  50-50  astrom_param_solved
c
c---------------------------------
 
      subroutine rdpmugbv_b50(ra,ra_err,
     .                        de,de_err,
     .                        pi,pi_err,
     .                        pr,pr_err,
     .                        pd,pd_err,
     .                        mg,mg_err,
     .                        mb,mb_err,
     .                        mr,mr_err,
     .                        sid,
     .                        apars,
     .                        b_50)
      implicit none
 
c
c input
c
      real*8    ra, ra_err
      real*8    de, de_err
      real*8    pi, pi_err
      real*8    pr, pr_err
      real*8    pd, pd_err
      real*8    mg, mg_err
      real*8    mb, mb_err
      real*8    mr, mr_err
      integer*8 sid
      integer*2 apars
 
c
c output
c
      byte      b_50(50)
 
c
c local
c
      integer i
      real*8 lra_err, lde_err
      real*8 lpi_err, lpr_err, lpd_err
      real*8 lmg_err, lmb_err, lmr_err
 
      integer*8 i8a_ra, i8b_ra, i8u_ra
      integer*8 i8a_de, i8b_de, i8u_de
      integer*8 i8a_pi, i8b_pi, i8u_pi
      integer*8 i8a_pr, i8b_pr, i8u_pr
      integer*8 i8a_pd, i8b_pd, i8u_pd
      integer*8 i8a_mg, i8b_mg, i8u_mg
      integer*8 i8a_mb, i8b_mb, i8u_mb
      integer*8 i8a_mr, i8b_mr, i8u_mr
 
      integer*8 i8_common
      byte      b8_common(8)
 
      real*8 piu, pru, pdu
 
      equivalence(i8_common,b8_common)
 
c
c-----------------------------------------------------
c
      i8a_ra = ra/360.0*17592186044416d0 + 0.5
      i8b_ra = 0
      lra_err = 0
      if (ra_err.gt.0) lra_err = log10(ra_err*1000)
      i8b_ra = 500*lra_err + 0.5
      if (i8b_ra.gt.4095) i8b_ra = 4095
 
      i8u_ra    = i8a_ra*4096 + i8b_ra
      i8_common = i8u_ra
 
      b_50(01) = b8_common(01)
      b_50(02) = b8_common(02)
      b_50(03) = b8_common(03)
      b_50(04) = b8_common(04)
      b_50(05) = b8_common(05)
      b_50(06) = b8_common(06)
      b_50(07) = b8_common(07)
 
c
c-----------------------------------------------------
c
      i8a_de = (de+90)/180.0*17592186044416d0 + 0.5
      i8b_de = 0
      lde_err = 0
      if (de_err.gt.0) lde_err = log10(de_err*1000)
      i8b_de = 500*lde_err + 0.5
      if (i8b_de.gt.4095) i8b_de = 4095
 
      i8u_de    = i8a_de*4096 + i8b_de
      i8_common = i8u_de
      b_50(08) = b8_common(01)
      b_50(09) = b8_common(02)
      b_50(10) = b8_common(03)
      b_50(11) = b8_common(04)
      b_50(12) = b8_common(05)
      b_50(13) = b8_common(06)
      b_50(14) = b8_common(07)
c
c-----------------------------------------------------
c
      piu = pi
      if (piu.lt.-1100d0) piu = -1100d0  ! alpha centauri is 0.768", 768 mas
      if (piu.gt. 1100d0) piu =  1100d0
 
      i8a_pi = (piu*100000.0d0 + 134217728.0d0) + 0.5
      if (i8a_pi.lt.000000000d0) i8a_pi = 000000000d0
      if (i8a_pi.gt.268435456d0) i8a_pi = 268435456d0
      i8b_pi = 0
      lpi_err = 0
      if (pi_err.gt.0) lpi_err = log10(pi_err*1000)
      i8b_pi = 500*lpi_err + 0.5
      if (i8b_pi.gt.4095) i8b_pi = 4095
 
      i8u_pi    = i8a_pi*4096 + i8b_pi
      i8_common = i8u_pi
      b_50(15) = b8_common(01)
      b_50(16) = b8_common(02)
      b_50(17) = b8_common(03)
      b_50(18) = b8_common(04)
      b_50(19) = b8_common(05)
c
c-----------------------------------------------------
c
      pru = pr
      if (pru.lt.-11000d0) pru = -11000d0  ! barnard's star is -10.4"/yr
      if (pru.gt. 11000d0) pru =  11000d0
 
      i8a_pr = 134217728 + pru*10000 + 0.5
      i8b_pr = 0
      lpr_err = 0
      if (pr_err.gt.0) lpr_err = log10(pr_err*1000)
      i8b_pr = 500*lpr_err + 0.5
      if (i8b_pr.gt.4095) i8b_pr = 4095
 
      i8u_pr    = i8a_pr*4096 + i8b_pr
      i8_common = i8u_pr
      b_50(20) = b8_common(01)
      b_50(21) = b8_common(02)
      b_50(22) = b8_common(03)
      b_50(23) = b8_common(04)
      b_50(24) = b8_common(05)
c
c-----------------------------------------------------
c
      pdu = pd
      if (pdu.lt.-11 000) pdu = -11000
      if (pdu.gt. 11 000) pdu =  11000
 
      i8a_pd = 134217728 + pdu*10000 + 0.5
      i8b_pd = 0
      lpd_err = 0
      if (pd_err.gt.0) lpd_err = log10(pd_err*1000)
      i8b_pd = 500*lpd_err + 0.5
      if (i8b_pd.gt.4095) i8b_pd = 4095
 
      i8u_pd    = i8a_pd*4096 + i8b_pd
      i8_common = i8u_pd
      b_50(25) = b8_common(01)
      b_50(26) = b8_common(02)
      b_50(27) = b8_common(03)
      b_50(28) = b8_common(04)
      b_50(29) = b8_common(05)
c
c-----------------------------------------------------
c
      i8a_mg = mg*40000 + 0.5
      i8b_mg = 0
      lmg_err = 0
      if (mg_err.gt.0) lmg_err = log10(mg_err*1000000)
      i8b_mg = 500*lmg_err + 0.5
      if (i8b_mg.gt.4095) i8b_mg = 4095
 
      i8u_mg    = i8a_mg*4096 + i8b_mg
      i8_common = i8u_mg
      b_50(30) = b8_common(01)
      b_50(31) = b8_common(02)
      b_50(32) = b8_common(03)
      b_50(33) = b8_common(04)
c
c-----------------------------------------------------
c
      i8a_mb = mb*40000 + 0.5
      i8b_mb = 0
      lmb_err = 0
      if (mb_err.gt.0) lmb_err = log10(mb_err*1000000)
      i8b_mb = 500*lmb_err + 0.5
      if (i8b_mb.gt.4095) i8b_mb = 4095
 
      i8u_mb    = i8a_mb*4096 + i8b_mb
      i8_common = i8u_mb
      b_50(34) = b8_common(01)
      b_50(35) = b8_common(02)
      b_50(36) = b8_common(03)
      b_50(37) = b8_common(04)
c
c-----------------------------------------------------
c
      i8a_mr = mr*40000 + 0.5
      i8b_mr = 0
      lmr_err = 0
      if (mr_err.gt.0) lmr_err = log10(mr_err*1000000)
      i8b_mr = 500*lmr_err + 0.5
      if (i8b_mr.gt.4095) i8b_mr = 4095
 
      i8u_mr    = i8a_mr*4096 + i8b_mr
      i8_common = i8u_mr
      b_50(38) = b8_common(01)
      b_50(39) = b8_common(02)
      b_50(40) = b8_common(03)
      b_50(41) = b8_common(04)
c
c-----------------------------------------------------
c
      i8_common = sid
      b_50(42) = b8_common(01)
      b_50(43) = b8_common(02)
      b_50(44) = b8_common(03)
      b_50(45) = b8_common(04)
      b_50(46) = b8_common(05)
      b_50(47) = b8_common(06)
      b_50(48) = b8_common(07)
      b_50(49) = b8_common(08)
c
c-----------------------------------------------------
c
      i8_common = apars
      b_50(50)  = b8_common(01)
c
c-----------------------------------------------------
c
 
 
c     write(44,144) ra, ra_err, lra_err, i8a_ra, i8b_ra,
c    .              de, de_err, lde_err, i8a_de, i8b_de
c     write(45,144) pi, pi_err, lpi_err, i8a_pi, i8b_pi,
c    .              pr, pr_err, lpr_err, i8a_pr, i8b_pr,
c    .              pd, pd_err, lpd_err, i8a_pd, i8b_pd
c     write(46,144) mg, mg_err, lmg_err, i8a_mg, i8b_mg,
c    .              mb, mb_err, lmb_err, i8a_mb, i8b_mb,
c    .              mr, mr_err, lmr_err, i8a_mr, i8b_mr
c     write(47,'(50Z2.2)') (b_50(i),i=1,50)
 
  144 format(1x,f14.10,f9.5,1x,f7.5,1x,i17,1x,i4,
     .       5x,f14.10,f9.5,1x,f7.5,1x,i17,1x,i4,
     .       5x,f14.10,f9.5,1x,f7.5,1x,i17,1x,i4)
 
      return
      end
 
 
 
 
c----------------------------------------------
c
c
c
c
      subroutine b50_rdpmugbv(b_50,
     .                        ra,ra_err,
     .                        de,de_err,
     .                        pi,pi_err,
     .                        pr,pr_err,
     .                        pd,pd_err,
     .                        mg,mg_err,
     .                        mb,mb_err,
     .                        mr,mr_err,
     .                        sid,
     .                        apars)
      implicit none
 
c
c intput
c
      byte      b_50(50)
 
c
c output
c
      real*8    ra, ra_err
      real*8    de, de_err
      real*8    pi, pi_err
      real*8    pr, pr_err
      real*8    pd, pd_err
      real*8    mg, mg_err
      real*8    mb, mb_err
      real*8    mr, mr_err
      integer*8 sid
      integer*2 apars
 
 
c
c local
c
      integer i
      real*8 lra_err, lde_err
      real*8 lpi_err, lpr_err, lpd_err
      real*8 lmg_err, lmb_err, lmr_err
 
      integer*8 i8a_ra, i8b_ra, i8u_ra
      integer*8 i8a_de, i8b_de, i8u_de
      integer*8 i8a_pi, i8b_pi, i8u_pi
      integer*8 i8a_pr, i8b_pr, i8u_pr
      integer*8 i8a_pd, i8b_pd, i8u_pd
      integer*8 i8a_mg, i8b_mg, i8u_mg
      integer*8 i8a_mb, i8b_mb, i8u_mb
      integer*8 i8a_mr, i8b_mr, i8u_mr
 
      integer*8 i8_common
      byte      b8_common(8)
 
 
      equivalence(i8_common,b8_common)
 
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(01)
      b8_common(02) = b_50(02)
      b8_common(03) = b_50(03)
      b8_common(04) = b_50(04)
      b8_common(05) = b_50(05)
      b8_common(06) = b_50(06)
      b8_common(07) = b_50(07)
      b8_common(08) = 0
 
      i8u_ra = i8_common
 
      i8a_ra = i8u_ra/4096
      i8b_ra = i8u_ra - i8u_ra/4096*4096
 
      ra = i8a_ra*360.0d0/17592186044416d0
 
      ra_err = 0.0
      if (i8b_ra.gt.0) then
          lra_err = i8b_ra/500.0
          ra_err  = 10**(lra_err)/1000.0
          endif
 
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(08)
      b8_common(02) = b_50(09)
      b8_common(03) = b_50(10)
      b8_common(04) = b_50(11)
      b8_common(05) = b_50(12)
      b8_common(06) = b_50(13)
      b8_common(07) = b_50(14)
      b8_common(08) = 0
 
      i8u_de = i8_common
 
      i8a_de = i8u_de/4096
      de = -90.0d0 + i8a_de*180/17592186044416d0
 
      i8b_de = i8u_de - i8u_de/4096*4096
      de_err = 0.0
      if (i8b_de.gt.0) then
          lde_err = i8b_de/500.0
          de_err  = 10**(lde_err)/1000.0
          endif
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(15)
      b8_common(02) = b_50(16)
      b8_common(03) = b_50(17)
      b8_common(04) = b_50(18)
      b8_common(05) = b_50(19)
      b8_common(06) = 0
      b8_common(07) = 0
      b8_common(08) = 0
 
      i8u_pi = i8_common
 
      i8a_pi = i8u_pi/4096
c     pi     = i8a_pi/100000.0d0
      pi     = (i8a_pi-134217728.0d0)/100000.0d0
 
      i8b_pi = i8u_pi - i8u_pi/4096*4096
      pi_err = 0.
      if (i8b_pi.gt.0) then
          lpi_err = i8b_pi/500.
          pi_err  = 10**(lpi_err)/1000.
          endif
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(20)
      b8_common(02) = b_50(21)
      b8_common(03) = b_50(22)
      b8_common(04) = b_50(23)
      b8_common(05) = b_50(24)
      b8_common(06) = 0
      b8_common(07) = 0
      b8_common(08) = 0
 
      i8u_pr = i8_common
      i8a_pr = i8u_pr/4096
      pr     = (i8a_pr-134217728)/10000.0
 
      i8b_pr = i8u_pr - i8u_pr/4096*4096
      pr_err  = 0.
      if (i8b_pr.gt.0) then
          lpr_err = i8b_pr/500.
          pr_err  = 10**(lpr_err)/1000.
          endif
c
c-----------------------------------------------------
c
 
      b8_common(01) = b_50(25)
      b8_common(02) = b_50(26)
      b8_common(03) = b_50(27)
      b8_common(04) = b_50(28)
      b8_common(05) = b_50(29)
      b8_common(06) = 0
      b8_common(07) = 0
      b8_common(08) = 0
 
      i8u_pd = i8_common
      i8a_pd = i8u_pd/4096
      pd     = (i8a_pd-134217728)/10000.0
 
      i8b_pd = i8u_pd - i8u_pd/4096*4096
      pd_err = 0
      if (i8b_pd.gt.0) then
          lpd_err = i8b_pd/500.
          pd_err  = 10**(lpd_err)/1000.
          endif
 
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(30)
      b8_common(02) = b_50(31)
      b8_common(03) = b_50(32)
      b8_common(04) = b_50(33)
      b8_common(05) = 0
      b8_common(06) = 0
      b8_common(07) = 0
      b8_common(08) = 0
 
      i8u_mg = i8_common
 
      i8a_mg = i8u_mg/4096
      mg     = i8a_mg/40000.0d0
 
      i8b_mg = i8u_mg - i8u_mg/4096*4096
      mg_err = 0.
      if (i8b_mg.gt.0) then
          lmg_err = i8b_mg/500.0
          mg_err  = 10**(lmg_err)/1000000d0
          endif
 
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(34)
      b8_common(02) = b_50(35)
      b8_common(03) = b_50(36)
      b8_common(04) = b_50(37)
      b8_common(05) = 0
      b8_common(06) = 0
      b8_common(07) = 0
      b8_common(08) = 0
 
      i8u_mb = i8_common
 
      i8a_mb = i8u_mb/4096
      mb     = i8a_mb/40000.0d0
 
      i8b_mb = i8u_mb - i8u_mb/4096*4096
      mb_err = 0.
      if (i8b_mb.gt.0) then
          lmb_err = i8b_mb/500.0
          mb_err  = 10**(lmb_err)/1000000d0
          endif
 
c
c-----------------------------------------------------
c
      b8_common(01) = b_50(38)
      b8_common(02) = b_50(39)
      b8_common(03) = b_50(40)
      b8_common(04) = b_50(41)
      b8_common(05) = 0
      b8_common(06) = 0
      b8_common(07) = 0
      b8_common(08) = 0
 
      i8u_mr = i8_common
 
      i8a_mr = i8u_mr/4096
      mr     = i8a_mr/40000.0d0
 
      i8b_mr = i8u_mr - i8u_mr/4096*4096
      mr_err = 0.
      if (i8b_mr.gt.0) then
          lmr_err = i8b_mr/500.0
          mr_err  = 10**(lmr_err)/1000000d0
          endif
 
c
c-----------------------------------------------------
c
 
      b8_common(01) = b_50(42)
      b8_common(02) = b_50(43)
      b8_common(03) = b_50(44)
      b8_common(04) = b_50(45)
      b8_common(05) = b_50(46)
      b8_common(06) = b_50(47)
      b8_common(07) = b_50(48)
      b8_common(08) = b_50(49)
      sid = i8_common
c
c-----------------------------------------------------
c
      apars = b_50(50)
 
      return
      end
 
 
c-----------------------------------------------------
c
c pack aux info#1
c
c
      subroutine rcor10_b10(ra_dec_corr,
     .                      ra_parallax_corr,
     .                      ra_pmra_corr,
     .                      ra_pmdec_corr,
     .                      dec_parallax_corr,
     .                      dec_pmra_corr,
     .                      dec_pmdec_corr,
     .                      parallax_pmra_corr,
     .                      parallax_pmdec_corr,
     .                      pmra_pmdec_corr,
     .                      b_10)
      implicit none
 
      real ra_dec_corr
      real ra_parallax_corr
      real ra_pmra_corr
      real ra_pmdec_corr
      real dec_parallax_corr
      real dec_pmra_corr
      real dec_pmdec_corr
      real parallax_pmra_corr
      real parallax_pmdec_corr
      real pmra_pmdec_corr
      byte b_10(10)
 
      if (abs(ra_dec_corr         ).gt.1.1) stop 'ra_de_corr bounds.'
      if (abs(ra_parallax_corr    ).gt.1.1) stop 'ra_pi_corr bounds.'
      if (abs(ra_pmra_corr        ).gt.1.1) stop 'ra_pr_corr bounds.'
      if (abs(ra_pmdec_corr       ).gt.1.1) stop 'ra_pd_corr bounds.'
      if (abs(dec_parallax_corr   ).gt.1.1) stop 'de_pi_corr bounds.'
      if (abs(dec_pmra_corr       ).gt.1.1) stop 'de_pr_corr bounds.'
      if (abs(dec_pmdec_corr      ).gt.1.1) stop 'de_pd_corr bounds.'
      if (abs(parallax_pmra_corr  ).gt.1.1) stop 'pi_pr_corr bounds.'
      if (abs(parallax_pmdec_corr ).gt.1.1) stop 'pi_pd_corr bounds.'
      if (abs(pmra_pmdec_corr     ).gt.1.1) stop 'pr_pd_corr bounds.'
 
      b_10(01) = int(100*ra_dec_corr        +1000+0.5) - 1000
      b_10(02) = int(100*ra_parallax_corr   +1000+0.5) - 1000
      b_10(03) = int(100*ra_pmra_corr       +1000+0.5) - 1000
      b_10(04) = int(100*ra_pmdec_corr      +1000+0.5) - 1000
      b_10(05) = int(100*dec_parallax_corr  +1000+0.5) - 1000
      b_10(06) = int(100*dec_pmra_corr      +1000+0.5) - 1000
      b_10(07) = int(100*dec_pmdec_corr     +1000+0.5) - 1000
      b_10(08) = int(100*parallax_pmra_corr +1000+0.5) - 1000
      b_10(09) = int(100*parallax_pmdec_corr+1000+0.5) - 1000
      b_10(10) = int(100*pmra_pmdec_corr    +1000+0.5) - 1000
 
      return
      end
 
 
 
c-----------------------------------------------------
c
c pack aux info#1
c
c
      subroutine b10_rcor10(b_10,
     .                      ra_dec_corr,
     .                      ra_parallax_corr,
     .                      ra_pmra_corr,
     .                      ra_pmdec_corr,
     .                      dec_parallax_corr,
     .                      dec_pmra_corr,
     .                      dec_pmdec_corr,
     .                      parallax_pmra_corr,
     .                      parallax_pmdec_corr,
     .                      pmra_pmdec_corr)
      implicit none
 
      byte b_10(10)
      real ra_dec_corr
      real ra_parallax_corr
      real ra_pmra_corr
      real ra_pmdec_corr
      real dec_parallax_corr
      real dec_pmra_corr
      real dec_pmdec_corr
      real parallax_pmra_corr
      real parallax_pmdec_corr
      real pmra_pmdec_corr
 
      ra_dec_corr         = b_10(01)/100.0
      ra_parallax_corr    = b_10(02)/100.0
      ra_pmra_corr        = b_10(03)/100.0
      ra_pmdec_corr       = b_10(04)/100.0
      dec_parallax_corr   = b_10(05)/100.0
      dec_pmra_corr       = b_10(06)/100.0
      dec_pmdec_corr      = b_10(07)/100.0
      parallax_pmra_corr  = b_10(08)/100.0
      parallax_pmdec_corr = b_10(09)/100.0
      pmra_pmdec_corr     = b_10(10)/100.0
 
      return
      end
 
 
c---------------------------------------------------------
c
c pack aux info#2
c
c
      subroutine astropar_b14(astrometric_n_obs_al,
     .                        astrometric_n_good_obs_al,
     .                        astrometric_n_bad_obs_al,
     .                        astrometric_gof_al,
     .                        astrometric_chi2_al,
     .                        astrometric_excess_noise,
     .                        astrometric_excess_noise_sig,
     .                        astrometric_primary_flag,
     .                        astrometric_matched_transits,
     .                        b_14)
 
      implicit none
 
      integer*2 astrometric_n_obs_al         !   1.3 bytes : uint    : 12-1564   : 11 bits : 2048
      integer*2 astrometric_n_good_obs_al    !   1.3 bytes : uint    : 10-1556   : 11 bits : 2048
      integer*2 astrometric_n_bad_obs_al     !   1.1 bytes : uint    :  0- 299   :  9 bits :  512
      real      astrometric_gof_al           !   2.0 bytes :  real*2 : -31-10382 : 16 bits :
      real      astrometric_chi2_al          !   2.0 bytes : ureal*2 : 2- 10000+ : 16 bits
      real      astrometric_excess_noise     !   2.0 bytes : ureal*2 : 0- 1032   : 16 bits
      real      astrometric_excess_noise_sig !   2.0 bytes : ureal82 : 0- 10000+ : 16 bits
c     integer*2 astrometric_params_solved    !   0.6 bytes : uint    : 2- 95     :  6 bits
      byte      astrometric_primary_flag     !                                      1 bit  :    2
      integer*2 astrometric_matched_transits !                       : 5-178     :  8 bits
      byte      b_14(14)
 
      byte   b1_2(2)
 
      integer*8 long
      byte      b8(8)
      equivalence (long,b8)
 
      long = astrometric_n_obs_al     *2048*512*2   ! 4 bytes
     .     + astrometric_n_good_obs_al     *512*2
     .     + astrometric_n_bad_obs_al          *2
     .     + astrometric_primary_flag
 
c     print*,'    1: ',astrometric_n_obs_al     *2048*512*2
c     print*,'    2: ',astrometric_n_good_obs_al     *512*2
c     print*,'    3: ',astrometric_n_bad_obs_al          *2
c     print*,'    4: ',astrometric_primary_flag
c     print*,' long: ',long
 
      b_14(01) = b8(01)
      b_14(02) = b8(02)
      b_14(03) = b8(03)
      b_14(04) = b8(04)
 
      call sr4_b2(astrometric_gof_al,b1_2)            ! signed 2-byte real
      b_14(05) = b1_2(1)
      b_14(06) = b1_2(2)
 
      call ur4_b2(astrometric_chi2_al,b1_2)
      b_14(07) = b1_2(1)
      b_14(08) = b1_2(2)
 
      call ur4_b2(astrometric_excess_noise,b1_2)
      b_14(09) = b1_2(1)
      b_14(10) = b1_2(2)
 
      call ur4_b2(astrometric_excess_noise_sig,b1_2)
      b_14(11) = b1_2(1)
      b_14(12) = b1_2(2)
 
      long = astrometric_matched_transits
      b_14(13) = b8(01)
      b_14(14) = b8(02)
 
      return
      end
 
c---------------------------------------------------------
c
c unpack aux info#2
c
c
      subroutine b14_astropar(b_14,
     .                        astrometric_n_obs_al,
     .                        astrometric_n_good_obs_al,
     .                        astrometric_n_bad_obs_al,
     .                        astrometric_gof_al,
     .                        astrometric_chi2_al,
     .                        astrometric_excess_noise,
     .                        astrometric_excess_noise_sig,
c    .                        astrometric_params_solved,      ! included in the basic line
     .                        astrometric_primary_flag,
     .                        astrometric_matched_transits)
 
      implicit none
 
      byte      b_14(14)
      integer*2 astrometric_n_obs_al         !   1.3 bytes : uint    : 12-1564   : 11 bits : 2048
      integer*2 astrometric_n_good_obs_al    !   1.3 bytes : uint    : 10-1556   : 11 bits : 2048
      integer*2 astrometric_n_bad_obs_al     !   1.1 bytes : uint    :  0- 299   :  9 bits :  512
      real      astrometric_gof_al           !   2.0 bytes :  real*2 : -31-10382 : 16 bits :
      real      astrometric_chi2_al          !   2.0 bytes : ureal*2 : 2- 10000+ : 16 bits
      real      astrometric_excess_noise     !   2.0 bytes : ureal*2 : 0- 1032   : 16 bits
      real      astrometric_excess_noise_sig !   2.0 bytes : ureal82 : 0- 10000+ : 16 bits
c     integer*2 astrometric_params_solved    !   0.6 bytes : uint    : 2- 95     :  6 bits
      byte      astrometric_primary_flag     !                                      1 bit  :    2
      integer*2 astrometric_matched_transits !                       : 5-178     :  8 bits
 
      byte   b1_2(2)
 
      integer*8 long
      byte      b8(8)
      equivalence (long,b8)
 
      long = 0
      b8(01) = b_14(01)
      b8(02) = b_14(02)
      b8(03) = b_14(03)
      b8(04) = b_14(04)
 
      astrometric_n_obs_al       = long / 2048/512/2
      astrometric_n_good_obs_al
     .        = (long-(long/2048/512/2)*2048*512*2)/(512*2)
      astrometric_n_bad_obs_al
     .        = (long-(long/512/2)*512*2)/2
      astrometric_primary_flag   = long-long/2*2
 
      b1_2(1) = b_14(05)
      b1_2(2) = b_14(06)
      call b2_sr4(b1_2,astrometric_gof_al)            ! signed 2-byte real
 
      b1_2(1) = b_14(07)
      b1_2(2) = b_14(08)
      call b2_ur4(b1_2,astrometric_chi2_al)
 
      b1_2(1) = b_14(09)
      b1_2(2) = b_14(10)
      call b2_ur4(b1_2,astrometric_excess_noise)
 
      b1_2(1) = b_14(11)
      b1_2(2) = b_14(12)
      call b2_ur4(b1_2,astrometric_excess_noise_sig)
 
      long = 0
      b8(01) = b_14(13)
      b8(02) = b_14(14)
      astrometric_matched_transits = long
 
 
      return
      end
 
 
c------------------------------------------------------------
c
c
c
      subroutine photpar_b11(phot_g_n_obs,
     .                       phot_bp_n_obs,
     .                       phot_rp_n_obs,
     .                       phot_bp_n_contaminated_transits,
     .                       phot_bp_n_blended_transits,
     .                       phot_rp_n_contaminated_transits,
     .                       phot_rp_n_blended_transits,
     .                       phot_bp_rp_excess_factor,
     .                       phot_proc_mode,
     .                       b_11)
 
      implicit none
 
      integer*2 phot_g_n_obs                    !  0 -  2395 : 12 bits (uint):  12 : 4096
      integer*2 phot_proc_mode                  !  0 -     2 :  2 bits (uint):   8
 
      integer*2 phot_bp_n_obs                   !  0 -   268 :  9 bits (uint):  21 :  512
      integer*2 phot_rp_n_obs                   !  0 -   263 :  9 bits (uint):  30 :  512
      integer*2 phot_bp_n_contaminated_transits !  0 -   211 :  8 bits (uint):  38 :  256
      integer*2 phot_bp_n_blended_transits      !  0 -   263 :  9 bits (uint):  47 :  512
      integer*2 phot_rp_n_contaminated_transits !  0 -   216 :  8 bits (uint):  55 :  256
      integer*2 phot_rp_n_blended_transits      !  0 -   262 :  9 bits (uint):  64 :  512
      real      phot_bp_rp_excess_factor        !  0 - 17495 :   2 bits (real): 16
 
      byte      b_11(11)
 
      integer*8 long
      byte      b8(8)
      equivalence (long,b8)
 
      integer*8 l_0512
      integer*8 l_4096
 
      byte b1_2(2)
 
      l_0512 =  512
      l_4096 = 4096
 
      long = phot_bp_n_obs          *l_0512*l_0512*l_0512*l_0512*l_0512
     .     + phot_rp_n_obs                 *l_0512*l_0512*l_0512*l_0512
     .     + phot_bp_n_contaminated_transits      *l_0512*l_0512*l_0512
     .     + phot_bp_n_blended_transits                  *l_0512*l_0512
     .     + phot_rp_n_contaminated_transits                    *l_0512
     .     + phot_rp_n_blended_transits
 
      b_11(01) = b8(01)
      b_11(02) = b8(02)
      b_11(03) = b8(03)
      b_11(04) = b8(04)
      b_11(05) = b8(05)
      b_11(06) = b8(06)
      b_11(07) = b8(07)
 
      long = phot_g_n_obs * 4 + phot_proc_mode
      b_11(08) = b8(01)
      b_11(09) = b8(02)
 
      call ur4_b2(phot_bp_rp_excess_factor,b1_2)
      b_11(10) = b1_2(01)
      b_11(11) = b1_2(02)
 
      return
      end
 
 
c------------------------------------------------------------
c
c
c
      subroutine b11_photpar(b_11,
     .                       phot_g_n_obs,
     .                       phot_bp_n_obs,
     .                       phot_rp_n_obs,
     .                       phot_bp_n_contaminated_transits,
     .                       phot_bp_n_blended_transits,
     .                       phot_rp_n_contaminated_transits,
     .                       phot_rp_n_blended_transits,
     .                       phot_bp_rp_excess_factor,
     .                       phot_proc_mode)
 
      implicit none
 
      byte      b_11(11)
 
      integer*2 phot_g_n_obs                    !  0 -  2395 : 12 bits (uint):  12 : 4096
      integer*2 phot_proc_mode                  !  0 -     2 :  2 bits (uint):   8
 
      integer*2 phot_bp_n_obs                   !  0 -   268 :  9 bits (uint):  21 :  512
      integer*2 phot_rp_n_obs                   !  0 -   263 :  9 bits (uint):  30 :  512
      integer*2 phot_bp_n_contaminated_transits !  0 -   211 :  8 bits (uint):  38 :  256
      integer*2 phot_bp_n_blended_transits      !  0 -   263 :  9 bits (uint):  47 :  512
      integer*2 phot_rp_n_contaminated_transits !  0 -   216 :  8 bits (uint):  55 :  256
      integer*2 phot_rp_n_blended_transits      !  0 -   262 :  9 bits (uint):  64 :  512
 
      real      phot_bp_rp_excess_factor        !  0 - 17495 :   2 bits (real): 16
 
      integer*8 long
      byte      b8(8)
      equivalence (long,b8)
 
      integer*8 l_0512
      integer*8 l_4096
 
      byte b1_2(2)
 
      l_0512 =  512
      l_4096 = 4096
 
      b8(01) = b_11(01)
      b8(02) = b_11(02)
      b8(03) = b_11(03)
      b8(04) = b_11(04)
      b8(05) = b_11(05)
      b8(06) = b_11(06)
      b8(07) = b_11(07)
      b8(08) = 0
 
 
      phot_bp_n_obs = long/l_0512**5
      phot_rp_n_obs = (long-(long/l_0512**5)*l_0512**5)/l_0512**4
      phot_bp_n_contaminated_transits
     .              = (long-(long/l_0512**4)*l_0512**4)/l_0512**3
      phot_bp_n_blended_transits
     .              = (long-(long/l_0512**3)*l_0512**3)/l_0512**2
      phot_rp_n_contaminated_transits
     .              = (long-(long/l_0512**2)*l_0512**2)/l_0512**1
      phot_rp_n_blended_transits
     .              = (long-(long/l_0512**1)*l_0512**1)/l_0512**0
 
      b8(01) = b_11(08)
      b8(02) = b_11(09)
      b8(03) = 0
      b8(04) = 0
      b8(05) = 0
      b8(06) = 0
      b8(07) = 0
      b8(08) = 0
      phot_g_n_obs   = long/4
      phot_proc_mode = long - long/4*4
 
      b1_2(01) = b_11(10)
      b1_2(02) = b_11(11)
      call b2_ur4(b1_2,phot_bp_rp_excess_factor)
 
 
      return
      end
 
 
 
 
c
c      program real2twobyte
c      implicit none
c
c      real*4  r4
c      real*4  r4o
c      real*4 ur4o
c      byte   b1_2(2)
c
c      integer q
c
c
c      do q = 01, 17
c
c         if (q.eq.01) r4 =       17495.0100
c         if (q.eq.02) r4 =        5235.4300
c         if (q.eq.03) r4 =           5.2433
c         if (q.eq.04) r4 =           0.9900
c         if (q.eq.05) r4 =           0.0888
c         if (q.eq.06) r4 =           0.0079
c         if (q.eq.07) r4 =      -17495.0100
c         if (q.eq.08) r4 =       -5235.4300
c         if (q.eq.09) r4 =          -5.2433
c         if (q.eq.10) r4 =          -0.9900
c         if (q.eq.11) r4 =          -0.0888
c         if (q.eq.12) r4 =          -0.0079
c         if (q.eq.13) r4 =   -17495010.0000
c         if (q.eq.14) r4 =   174950100.0000
c         if (q.eq.15) r4 =  -174950100.0000
c         if (q.eq.16) r4 =  1749501000.0000
c         if (q.eq.17) r4 = -1749501000.0000
c         call ur4_b2(r4,b1_2)
c         call b2_ur4(b1_2,ur4o)
c
c         call sr4_b2(r4,b1_2)
c         call b2_sr4(b1_2,r4o)
c
c
c         write(*,119) q,b1_2(1),b1_2(2),r4,r4o,ur4o
c  119    format(1x,i3,1x,i4,1x,i4,1x,f20.7,1x,f20.7,1x,f20.7)
c
c         enddo
c
c      stop
c      end
c
 
 
c------------------------------
c
c
      subroutine sr4_b2(r,b1_2)
      implicit none
 
      real r
      byte b1_2(2)
 
      real     lr
      integer ilr
      real      e
      real     mant
      integer imant
 
      integer i1
      integer i2
 
      real ar
 
      integer*2 short
      byte      bb(2)
      equivalence(short,bb)
 
 
      ar = abs(r)
 
      lr = log10(ar)
      ilr = int(lr)
 
      if (lr.lt.0) ilr = ilr -1
 
      i2 = ilr
      if (i2.gt. 7) i2 =  7
      if (i2.lt.-7) i2 = -7
 
      e = 10**(i2+1.0)
 
      mant = r/e
 
      if (mant.gt. 0.999) mant =  0.999
      if (mant.lt.-0.999) mant = -0.999
 
      imant = int(5000 + mant*2000 + 0.5) - 5000
 
      if (imant.le.-2047) stop 'imant.le.-2047'
      if (imant.ge. 2047) stop 'imant.le. 2047'
 
      i1 = imant
 
      short = i1*16 + (i2+7)
      if (i1.lt.0) short = i1*16 - (i2+7)
 
      b1_2(1) = bb(1)
      b1_2(2) = bb(2)
 
c     write(*,119) r, lr, ilr, e, mant, imant, i2,
c    .             (i1/200.0) * 10**(i2*1.0),
c    .             short, b1_2(1), b1_2(2)
  119 format(1x,f22.9,1x,f6.3,1x,i5,1x,f22.9,
     .       1x,f9.6,1x,i5,1x,i3,
     .       1x,f22.9,
     .       1x,i6,1x,i4,1x,i4)
 
      return
      end
 
 
 
c------------------------------
c
c
      subroutine b2_sr4(b1_2,r)
      implicit none
 
      byte b1_2(2)
      real r
 
      real     lr
      integer ilr
      real      e
      real     mant
      integer imant
 
      integer i1
      integer i2
 
      integer*2 short
      byte      bb(2)
      equivalence(short,bb)
 
      bb(1) = b1_2(1)
      bb(2) = b1_2(2)
 
      i1    = short/16
      i2    = short - i1*16 - 7
      if (i1.lt.0) i2 = -(short - i1*16) - 7
 
      mant = i1/2000.0d0
 
      e = 10**(i2+1.0)
 
      r = mant*e
 
      return
      end
 
c------------------------------
c
c
      subroutine ur4_b2(r,b1_2)
      implicit none
 
      real r
      byte b1_2(2)
 
      real     lr
      integer ilr
      real      e
      real     mant
      integer imant
 
      integer i1
      integer i2
 
      real ar
 
      integer*4 long
      byte      bbbb(4)
      equivalence(long,bbbb)
 
 
      ar = abs(r)
 
      lr = log10(ar)
      ilr = int(lr)
 
      if (lr.lt.0) ilr = ilr -1
 
      i2 = ilr
      if (i2.gt. 7) i2 =  7
      if (i2.lt.-7) i2 = -7
 
      e = 10**(i2+1.0)
 
      mant = r/e
 
      if (mant.gt. 0.999) mant =  0.999
      if (mant.lt. 0.000) mant =  0.000
 
      imant = int(mant*4000+0.5)
 
      if (imant.lt.    0) stop 'imant.lt.    0'
      if (imant.gt. 4090) stop 'imant.le. 4090'
 
      i1 = imant
 
      long = i1*16 + (i2+7)
 
      b1_2(1) = bbbb(1)
      b1_2(2) = bbbb(2)
 
c     write(*,119) r, lr, ilr, e, mant, imant, i2,
c    .             (i1/200.0) * 10**(i2*1.0),
c    .             short, b1_2(1), b1_2(2)
  119 format(1x,f22.9,1x,f6.3,1x,i5,1x,f22.9,
     .       1x,f9.6,1x,i5,1x,i3,
     .       1x,f22.9,
     .       1x,i6,1x,i4,1x,i4)
 
      return
      end
 
 
 
c------------------------------
c
c
      subroutine b2_ur4(b1_2,r)
      implicit none
 
      byte b1_2(2)
      real r
 
      real     lr
      integer ilr
      real      e
      real     mant
      integer imant
 
      integer i1
      integer i2
 
      integer*4 long
      byte      bbbb(4)
      equivalence(long,bbbb)
 
      bbbb(1) = b1_2(1)
      bbbb(2) = b1_2(2)
      bbbb(3) = 0
      bbbb(4) = 0
 
      i1    = long/16
      i2    = long - i1*16 - 7
 
      mant = i1/4000.0d0
 
      e = 10**(i2+1.0)
 
      r = mant*e
 
      return
      end
 
 
 
 
c------------------------------------
c
c
c
c
      subroutine showfull_b85(b_85,TF,UNIT)
      implicit none
 
      byte    b_85(85)
      logical TF
      integer UNIT
 
      real*8    ra, ra_err
      real*8    de, de_err
      real*8    pi, pi_err
      real*8    pr, pr_err
      real*8    pd, pd_err
      real*8    mg, mg_err
      real*8    mb, mb_err
      real*8    mr, mr_err
      integer*8 sid
      integer*2 apars
 
      real ra_dec_corr
      real ra_parallax_corr
      real ra_pmra_corr
      real ra_pmdec_corr
      real dec_parallax_corr
      real dec_pmra_corr
      real dec_pmdec_corr
      real parallax_pmra_corr
      real parallax_pmdec_corr
      real pmra_pmdec_corr
 
      byte      b_14(14)
      integer*2 astrometric_n_obs_al
      integer*2 astrometric_n_good_obs_al
      integer*2 astrometric_n_bad_obs_al
      real      astrometric_gof_al
      real      astrometric_chi2_al
      real      astrometric_excess_noise
      real      astrometric_excess_noise_sig
      byte      astrometric_primary_flag
      integer*2 astrometric_matched_transits
 
      integer*2 phot_g_n_obs
      integer*2 phot_proc_mode
      integer*2 phot_bp_n_obs
      integer*2 phot_rp_n_obs
      integer*2 phot_bp_n_contaminated_transits
      integer*2 phot_bp_n_blended_transits
      integer*2 phot_rp_n_contaminated_transits
      integer*2 phot_rp_n_blended_transits
      real      phot_bp_rp_excess_factor
      real      ruwe_id
 
      integer   isrv
 
      logical   FIRST
      data FIRST/.true./
      common /FIRSTo/FIRST
 
 
      call b50_rdpmugbv(b_85,
     .                  ra,ra_err,
     .                  de,de_err,
     .                  pi,pi_err,
     .                  pr,pr_err,
     .                  pd,pd_err,
     .                  mg,mg_err,
     .                  mb,mb_err,
     .                  mr,mr_err,
     .                  sid,
     .                  apars)
 
 
      call b10_rcor10(b_85(51),
     .                ra_dec_corr,
     .                ra_parallax_corr,
     .                ra_pmra_corr,
     .                ra_pmdec_corr,
     .                dec_parallax_corr,
     .                dec_pmra_corr,
     .                dec_pmdec_corr,
     .                parallax_pmra_corr,
     .                parallax_pmdec_corr,
     .                pmra_pmdec_corr)
 
 
      call b14_astropar(b_85(61),
     .                  astrometric_n_obs_al,
     .                  astrometric_n_good_obs_al,
     .                  astrometric_n_bad_obs_al,
     .                  astrometric_gof_al,
     .                  astrometric_chi2_al,
     .                  astrometric_excess_noise,
     .                  astrometric_excess_noise_sig,
     .                  astrometric_primary_flag,
     .                  astrometric_matched_transits)
 
 
      call b11_photpar(b_85(75),
     .                 phot_g_n_obs,
     .                 phot_bp_n_obs,
     .                 phot_rp_n_obs,
     .                 phot_bp_n_contaminated_transits,
     .                 phot_bp_n_blended_transits,
     .                 phot_rp_n_contaminated_transits,
     .                 phot_rp_n_blended_transits,
     .                 phot_bp_rp_excess_factor,
     .                 phot_proc_mode)
 
 
      isrv = 0
      if (apars.lt.0) then
          isrv = 1
          apars = -apars
          endif
 
c     write(51,121) ra, de,
c    .              ra_err, de_err,
c    .              mg, mg_err,
c    .              mb, mb_err,
c    .              mr, mr_err,
c    .              pr, pr_err,
c    .              pd, pd_err,
c    .              pi, pi_err,
c    .              isrv,
c    .              sid
  121 format(f14.10,1x,f14.10,1x,
     .       f7.3,1x,f7.3,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f8.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       i1,1x,i19)
  115 format('#    (deg)    ',1x,'      (deg)   ',1x,
     .       '  (mas)',1x,'  (mas)',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       ' (mag)',1x,' (mag) ',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '(mas/yr)',1x,'(m/y)',1x,
     .       '   (mas)',1x,'(mas)',1x,
     .       '                    ')
  116 format('#RIGHT ASCENSN',1x,'   DECLINATION',1x,
     .       '   eRA ',1x,'   eDEC',1x,
     .       '  GMAG',1x,' eGMAG ',1x,
     .       '  BMAG',1x,' eBMAG ',1x,
     .       '  RMAG',1x,' eRMAG ',1x,
     .       '    muRA',1x,'emuRA',1x,
     .       '   muDEC',1x,'emuDE',1x,
     .       '     PI ',1x,'  ePI',1x,
     .       ' GAIA_DR2_ID_NUMBER ')
  117 format('#.............',1x,'..............',1x,
     .       '.......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '......',1x,'.......',1x,
     .       '......',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '........',1x,'.....',1x,
     .       '....................')
c
c----------------------------------------------------
c
c     write(52,123) ra_dec_corr,
c    .              ra_parallax_corr,
c    .              ra_pmra_corr,
c    .              ra_pmdec_corr,
c    .              dec_parallax_corr,
c    .              dec_pmra_corr,
c    .              dec_pmdec_corr,
c    .              parallax_pmra_corr,
c    .              parallax_pmdec_corr,
c    .              pmra_pmdec_corr,
c    .              sid
  123 format('   ',10(f6.2,1x),i19)
  423 format('#  ',' ra_dec   ',
     .             ' ra_par   ',
     .             ' ra_pmra  ',
     .             ' ra_pmde  ',
     .             ' de_par   ',
     .             ' de_pmra  ',
     .             ' de_pmde  ',
     .             ' par_pmra ',
     .             ' par_pmde ',
     .             ' star_ID    ')
c
c----------------------------------------------------
c
c     write(53,124) astrometric_n_obs_al,
c    .              astrometric_n_good_obs_al,
c    .              astrometric_n_bad_obs_al,
c    .              astrometric_gof_al,
c    .              astrometric_chi2_al,
c    .              astrometric_excess_noise,
c    .              astrometric_excess_noise_sig,
c    .              apars,
c    .              astrometric_primary_flag,
c    .              astrometric_matched_transits,
c    .              sid
  124 format('   ',i4,1x,i4,1x,i4,
     .          1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,
     .          1x,i4,1x,i4,1x,i4,1x,i19)
  424 format('#  ','n_ob',1x,'n_go',1x,'n_ba',
     .                    1x,'    gof_al',
     .                    1x,'   chi2_al',
     .                    1x,'    xs_noi',
     .                    1x,' xs_noisig',
     .                    1x,'psol',1x,'prim',1x,'matt',
     .                    1x,' dr3_star_id ')
c
c----------------------------------------------------
c
c     write(54,125) phot_g_n_obs,
c    .              phot_bp_n_obs,
c    .              phot_rp_n_obs,
c    .              phot_bp_n_contaminated_transits,
c    .              phot_bp_n_blended_transits,
c    .              phot_rp_n_contaminated_transits,
c    .              phot_rp_n_blended_transits,
c    .              phot_proc_mode,
c    .              phot_bp_rp_excess_factor,
c    .              sid
  125 format('   ',i4,1x,i4,1x,i4,
     .       1x,i4,1x,i4,1x,i4,1x,i4,
     .       1x,i4,1x,f8.2,1x,i19)
 
      if (FIRST) then
         write(UNIT,504)
         write(UNIT,502)
         write(UNIT,501)
         write(UNIT,503)
         write(UNIT,502)
         FIRST = .false.
         endif
      write(UNIT,500) ra, de,
     .              ra_err, de_err,
     .              mg, mg_err,
     .              mb, mb_err,
     .              mr, mr_err,
     .              pr, pr_err,
     .              pd, pd_err,
     .              pi, pi_err,
     .              isrv,
     .              sid,
     .              ra_dec_corr,
     .              ra_parallax_corr,
     .              ra_pmra_corr,
     .              ra_pmdec_corr,
     .              dec_parallax_corr,
     .              dec_pmra_corr,
     .              dec_pmdec_corr,
     .              parallax_pmra_corr,
     .              parallax_pmdec_corr,
     .              pmra_pmdec_corr,
     .              astrometric_n_obs_al,
     .              astrometric_n_good_obs_al,
     .              astrometric_n_bad_obs_al,
     .              astrometric_gof_al,
     .              astrometric_chi2_al,
     .              astrometric_excess_noise,
     .              astrometric_excess_noise_sig,
     .              apars,
     .              astrometric_primary_flag,
     .              astrometric_matched_transits,
     .              phot_g_n_obs,
     .              phot_bp_n_obs,
     .              phot_rp_n_obs,
     .              phot_bp_n_contaminated_transits,
     .              phot_bp_n_blended_transits,
     .              phot_rp_n_contaminated_transits,
     .              phot_rp_n_blended_transits,
     .              phot_proc_mode,
     .              phot_bp_rp_excess_factor,
     .              ruwe_id(sid)
 
  500 format(f14.10,1x,f14.10,1x,
     .       f7.3,1x,f7.3,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f6.3,1x,f7.5,1x,
     .       f8.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       f8.3,1x,f5.3,1x,
     .       i1,1x,i19,
     .       1x,f6.2,9(1x,f6.2),
     .       1x,i4,1x,i4,1x,i4,
     .       1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,
     .       1x,i4,1x,i4,1x,i4,
     .       1x,i4,1x,i4,1x,i4,
     .       1x,i4,1x,i4,1x,i4,1x,i4,
     .       1x,i4,1x,f8.2,1x,f6.2)
 
  501 format('#   RA (deg)  ',1x,'   DEC (deg)  ',1x,
     .       1x,'   eRA ',1x,'  eDEC',
     .       1x,'  GMAG',1x,' eGMAG ',
     .       1x,'  BMAG',1x,' eBMAG ',
     .       1x,'  RMAG',1x,' eRMAG ',
     .       1x,'    muRA',1x,'emuRA',
     .       1x,'   muDEC',1x,'emuDE',
     .       1x,'     PI ',1x,' ePI',
     .       1x,'RV?',
     .       1x,'GAIA_DR2_ID_NUMBER',
     .       1x,'c_rade',
     .       1x,'c_rapi',
     .       1x,'c_rapr',
     .       1x,'c_rapd',
     .       1x,'c_depi',
     .       1x,'c_depr',
     .       1x,'c_depd',
     .       1x,'c_pipr',
     .       1x,'c_pipd',
     .       1x,'c_prpd',
     .       1x,'n_ob',
     .       1x,'n_go',
     .       1x,'n_ba',
     .       1x,'    gof_al',
     .       1x,'   chi2_al',
     .       1x,'    xs_noi',
     .       1x,' xs_noisig',
     .       1x,'psol',1x,'prim',1x,'matt',
     .       1x,'g_no',1x,'b_no',1x,'r_no',
     .       1x,'b_co',1x,'b_bo',
     .       1x,'r_co',1x,'r_bo',
     .       1x,'proc',1x,'bprp_xss')
 
  502 format('#-------------',1x,'--------------',1x,
     .       1x,'-------',1x,'------',
     .       1x,'------',1x,'-------',
     .       1x,'------',1x,'-------',
     .       1x,'------',1x,'-------',
     .       1x,'--------',1x,'-----',
     .       1x,'--------',1x,'-----',
     .       1x,'--------',1x,'-----',
     .       1x,'-',
     .       1x,'-------------------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'------',
     .       1x,'----',
     .       1x,'----',
     .       1x,'----',
     .       1x,'----------',
     .       1x,'----------',
     .       1x,'----------',
     .       1x,'----------',
     .       1x,'----',1x,'----',1x,'----',
     .       1x,'----',1x,'----',1x,'----',
     .       1x,'----',1x,'----',
     .       1x,'----',1x,'----',
     .       1x,'----',1x,'--------')
 
  503 format('#    (01)     ',1x,'    (02)      ',1x,
     .       1x,'  (03) ',1x,' (04) ',
     .       1x,' (05) ',1x,'  (06) ',
     .       1x,' (07) ',1x,'  (08) ',
     .       1x,' (09) ',1x,'  (10) ',
     .       1x,'   (11) ',1x,' (12)',
     .       1x,'   (13) ',1x,' (14)',
     .       1x,'   (15) ',1x,' (16)',
     .       1x,' ',
     .       1x,'     (18)          ',
     .       1x,'  (19)',
     .       1x,'  (20)',
     .       1x,'  (21)',
     .       1x,'  (22)',
     .       1x,'  (23)',
     .       1x,'  (24)',
     .       1x,'  (25)',
     .       1x,'  (26)',
     .       1x,'  (27)',
     .       1x,'  (28)',
     .       1x,'(29)',
     .       1x,'(30)',
     .       1x,'(31)',
     .       1x,'     (32) ',
     .       1x,'     (33) ',
     .       1x,'     (34) ',
     .       1x,'     (35) ',
     .       1x,'(36)',1x,'(37)',1x,'(38)',
     .       1x,'(39)',1x,'(40)',1x,'(41)',
     .       1x,'(42)',1x,'(43)',
     .       1x,'(44)',1x,'(45)',
     .       1x,'(46)',1x,'   (47) ')
 
  504 format('#                                                ',/,
     .       '# (01) RA         RIGHT ASCENSION      -- degrees',/,
     .       '# (02) DEC        DECLINATION          -- degrees',/,
     .       '# (03) eRA        RA ERROR             -- mas    ',/,
     .       '# (04) eDE        DEC ERROR            -- mas    ',/,
     .       '# (05) GMAG       GMAG                 -- mags   ',/,
     .       '# (06) eGMAG      GMAG ERROR           -- mags   ',/,
     .       '# (07) BMAG       BMAG                 -- mags   ',/,
     .       '# (08) eBMAG      BMAG ERROR           -- mags   ',/,
     .       '# (09) RMAG       RMAG                 -- mags   ',/,
     .       '# (10) eRMAG      RMAG ERROR           -- mags   ',/,
     .       '# (11) PMRA       PM RA                -- mas/yr ',/,
     .       '# (12) ePMRA      ERROR PM RA          -- mas/yr ',/,
     .       '# (13) PMDE       PM DEC               -- mas/yr ',/,
     .       '# (14) ePMDE      ERROR PM DEC         -- mas/yr ',/,
     .       '# (15) PI         PARALLAX             -- mas    ',/,
     .       '# (16) ePI        ERROR PARALLAX       -- mas    ',/,
     .       '# (17) rv?        IS THERE AN RV?      -- 0 or 1 ',/,
     .       '# (18) ID         GAIA DR3 SOURCE ID   -- long  ',/,
     .       '# (19) c_rade     CORREL RA-DEC        -- -1 to 1 ',/,
     .       '# (20) c_rapi     CORREL RA-PARALLAX   -- -1 to 1 ',/,
     .       '# (21) c_rapr     CORREL RA-PMRA       -- -1 to 1 ',/,
     .       '# (22) c_rapd     CORREL RA-PMDEC      -- -1 to 1 ',/,
     .       '# (23) c_depi     CORREL DEC-PARALLAX  -- -1 to 1 ',/,
     .       '# (24) c_depr     CORREL DEC-PMRA      -- -1 to 1 ',/,
     .       '# (25) c_depd     CORREL DEC-PMDEC     -- -1 to 1 ',/,
     .       '# (26) c_pipr     CORREL PARALLAX-PMRA -- -1 to 1 ',/,
     .       '# (27) c_pipd     CORREL PARALLAX-PMDE -- -1 to 1 ',/,
     .       '# (28) c_prpd     CORREL PMRA-PMDE     -- -1 to 1 ',/,
     .       '# (29) n_ob       NUMBER OF OBSERVATNs -- short   ',/,
     .       '# (30) n_go       NUMBER OF GOOD OBSNS -- short   ',/,
     .       '# (31) n_ba       NUMBER OF BAD OBSNS  -- shor    ',/,
     .       '# (32) gof_al     ASTROMETRIC GOF      --         ',/,
     .       '# (33) chi2_al    ASTROMETRIC CHIS2    --         ',/,
     .       '# (34) xs_noi     EXCESS NOISE         --         ',/,
     .       '# (35) xs_noisig  EXCESS NOISE SIG     --         ',/,
     .       '# (36) psol       ASTROM PARAM SOLVED  -- short   ',/,
     .       '# (37) prim       PRIMARY FLAG         -- 1 or 0  ',/,
     .       '# (38) matt       ASTROM MATCHD TRANSs -- short   ',/,
     .       '# (39) g_no       NUM OBSs IN GBAND    -- short   ',/,
     .       '# (40) b_no       NUM OBSs IN BBAND    -- short   ',/,
     .       '# (41) r_no       NUM OBSs IN RBAND    -- short   ',/,
     .       '# (42) b_co       NUM CONTAM IN BBAND  -- short   ',/,
     .       '# (43) b_bo       NUM BLENDED IN BBAND -- short   ',/,
     .       '# (44) r_co       NUM CONTAM IN RBAND  -- short   ',/,
     .       '# (45) r_bo       NUM BLENDED IN RBAND -- short   ',/,
     .       '# (46) proc       PROC PHOT MODE 0=GOLD - 0, 1, 2 ',/,
     .       '# (47) bprp_xss   EXCESS BP/RP                    ',/,
     .       '#                                                 ')
 
      if (TF) then
         write(UNIT,502)
         write(UNIT,503)
         write(UNIT,501)
         write(UNIT,502)
         FIRST = .false.
         endif
 
      return
      end
 
 
c---------------------------------------------------------------
c
c FROM: ../GAIA_DR3_ASSEMBLE/11.GaiaSelect_RUWE/prog7_function.F
c
c---------------------------------------------------------------
 
 
c
c      program return_ruwe
c      implicit none
c
c      real ruwe_id
c      integer*8 ID
c      real ruwe
c
c      integer N
c      character*21 STRING20
c
c      N = 0
c    1 continue
c      read(*,*,end=2) STRING20
c      read(STRING20,*,end=2) ID
c      N = N + 1
cc     write(*,'(i6,1x,a20,1x,i20,1x,f9.3)') N, STRING20, ID
c      ruwe = ruwe_id(ID)
c      write(*,'(i6,1x,a20,1x,i20,1x,f9.3)') N, STRING20, ID, ruwe
c      goto 1
c    2 continue
c
c      stop
c      end
c
c
 
 
      real function ruwe_id(IDi)
      implicit none
 
      integer*8 IDi
 
      byte b_10(10)
      byte b_08(08)
      byte b_02(02)
      integer*8 i8, i8o
      integer*2 i2
 
      equivalence(i8,b_08)
      equivalence(i2,b_02)
 
      byte b_10_R(10,1000)
 
      integer*4 Os_4
      integer*8 Os_8
      integer*8 Os_8x10
 
      integer   i
      integer   L, Ls
 
      integer*8 R, Rs, Bs
      integer*8 ID
 
      logical   first
      integer*8 id_r(1811710)
      common / common_IDR1 / first
      common / common_IDR2 / id_r
      data first / .true. /
 
      integer   NIT
      integer*8 R1, R2, RR
 
      integer*8 L1, L2, LL
 
      real      ruwe
 
 
      i8 = 0
      i2 = 0
      Rs = 1811710
 
      if (first) then
          open(98,file=
     .   'NONE',
     .            status='old',
     .            err=900,
     .            recl=14493680,
     .            form='UNFORMATTED',
     .            access='DIRECT')
          read(98,rec=1251) id_r
          close(98)
          first = .false.
          if (.false.) then
              print*,' '
              print*,'FIRST! '
              print*,' '
              write(*,'(2i8,1x,i20)') 0, 1, id_r(000001)
              do i = 1, 10
                 R = 1811710*i/10
                 write(*,'(2i8,1x,i20)') i, R, id_r(R)
                 enddo
              print*,' '
              endif
          endif
 
 
 
 
      ID = IDi
 
      if (id.ge.id_r(1811710)) then
          R = 1811710
          goto 13
          endif
 
      r1 = 0000001
      r2 = 1811709
 
      NIT = 0
   11 continue
      NIT = NIT + 1
      rr = (r1+r2)/2
c     write(*,555) NIT,R1,RR,R2,
c    .             id_r(R1  ), id_r(RR),
c    .             id_r(RR+1), id_r(R2), ID
  555 format(1x,i3,5x,i7,1x,i7,1x,i7,5x,
     .                i20,1x,i20,1x,
     .                i20,1x,i20,5x,i20)
      if (NIT.gt.25) stop
      if (ID.lt.id_r(RR)) then
          r2 = rr
          goto 11
          endif
      if (ID.gt.id_r(RR+1)) then
          r1 = rr
          goto 11
          endif
   13 continue
c     print*,' '
c     write(*,'('' R1  : '',i7,1x,i20)') R1  ,id_r(R1)
c     write(*,'('' RR  : '',i7,1x,i20)') RR  ,id_r(RR)
c     write(*,'(''     : '',i7,1x,i20)')  0  ,ID
c     write(*,'('' RR+1: '',i7,1x,i20)') RR+1,id_r(RR+1)
c     write(*,'('' R2  : '',i7,1x,i20)') R2  ,id_r(R2)
c     print*,' '
 
      open(98,file=
     .   'NONE',
     .        status='old',
     .        err=900,
     .        recl=10000,
     .        form='UNFORMATTED',
     .        access='DIRECT')
 
      read(98,rec=RR) b_10_R
 
      ruwe = -9.99
 
      LL = 0
      do L = 1, 1000
         b_08(01) = b_10_R(01,L)
         b_08(02) = b_10_R(02,L)
         b_08(03) = b_10_R(03,L)
         b_08(04) = b_10_R(04,L)
         b_08(05) = b_10_R(05,L)
         b_08(06) = b_10_R(06,L)
         b_08(07) = b_10_R(07,L)
         b_08(08) = b_10_R(08,L)
         b_02(01) = b_10_R(09,L)
         b_02(02) = b_10_R(10,L)
         if (i8.eq.ID) then
             LL = L
c            write(*,'(i7,1x,i4,1x,i20,1x,i6,5x,10i5)')
c    .          RR,L,i8,i2,(b_10_R(i,L),i=1,10)
             ruwe = i2/100.0
             endif
         enddo
 
c     if (LL.eq.0) stop 'NOT FOUND'
 
      ruwe_id = ruwe
 
      return
  900 print*,'file open error...'
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/ROUTINES/TRANS/glob_fit6nrDPU.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine glob_fit6nrDPU(x1,y1,x2,y2,uu,Ns,
     .                          A,B,C,D,x1o,y1o,x2o,y2o)
      implicit none
 
      real*8 x1(*), y1(*)
      real*8 x2(*), y2(*)
      integer uu(*)
      integer Ns
      real*8 A, B, C, D
      real*8 x1o, y1o, x2o, y2o
 
      real*8 sxx, sx, swx, szx
      real*8 syy, sy, swy, szy
      real*8 sw,  sz, sxy
      real*8 dlta
      real*8 dsxx, dsyy, dsxy
      real*8 dswx, dswy, dszx, dszy
 
      integer n, nnn, NNu
 
      real*8 x1a, y1a, x1b, y1b
      real*8 x2a, y2a, x2b, y2b
      real*8 x2ae,y2ae
      real*8 dx1, dy1
      real*8 dx2, dy2
 
      if (Ns.lt.3) goto 999
 
      x1a = 0.
      x1b = 0.
      y1a = 0.
      y1b = 0.
 
      x2a = 0.
      x2b = 0.
      y2a = 0.
      y2b = 0.
 
      x1o = 0
      y1o = 0
      x2o = 0
      y2o = 0
      nnn = 0
      do n = 1, Ns
           x1o = x1o + uu(n)*x1(n)
           y1o = y1o + uu(n)*y1(n)
           x2o = x2o + uu(n)*x2(n)
           y2o = y2o + uu(n)*y2(n)
           nnn = nnn + uu(n)
           enddo
 
      x1o = x1o/NNN
      y1o = y1o/NNN
      x2o = x2o/NNN
      y2o = y2o/NNN
 
      if (NNN.eq.1) then
         A = 1.000
         B = 0.000
         C = 0.000
         D = 1.000
         return
         endif
 
      if (NNN.eq.2) then
          NNu = 0
          do N = 1, Ns
             if (uu(N).eq.1.and.NNu.eq.2) stop 'prob:  NNu.eq.2'
             if (uu(N).eq.1.and.NNu.eq.1) then
                 x1b = x1(n)
                 y1b = y1(n)
                 x2b = x2(n)
                 y2b = y2(n)
                 NNu = 2
                 endif
             if (uu(N).eq.1.and.NNu.eq.0) then
                 x1a = x1(n)
                 y1a = y1(n)
                 x2a = x2(n)
                 y2a = y2(n)
                 NNu = 1
                 endif
             enddo
          x1o = (x1a+x1b)/2.0
          y1o = (y1a+y1b)/2.0
          x2o = (x2a+x2b)/2.0
          y2o = (y2a+y2b)/2.0
          dx1 = (x1b-x1a)
          dy1 = (y1b-y1a)
          dx2 = (x2b-x2a)
          dy2 = (y2b-y2a)
          A =  (dx1*dx2+dy1*dy2)/(dx1**2+dy1**2)
          B =  (dy1*dx2-dx1*dy2)/(dx1**2+dy1**2)
          C = -B
          D =  A
          x2ae = x2o + A*(x1a-x1o) + B*(y1a-y1o)
          y2ae = y2o + C*(x1a-x1o) + D*(y1a-y1o)
          return
          endif
 
 
      sxx = 0.0
      sx  = 0.0
      syy = 0.0
      sy  = 0.0
      swx = 0.0
      swy = 0.0
      szx = 0.0
      szy = 0.0
      sw  = 0.0
      sz  = 0.0
      sxy = 0.0
      do n = 1, Ns
           sxy = sxy + uu(n)*(x1(n)-x1o)*(y1(n)-y1o)
           sxx = sxx + uu(n)*(x1(n)-x1o)*(x1(n)-x1o)
           sx  = sx  + uu(n)*(x1(n)-x1o)
           syy = syy + uu(n)*(y1(n)-y1o)*(y1(n)-y1o)
           sy  = sy  + uu(n)*(y1(n)-y1o)
           swx = swx + uu(n)*(x2(n)-x2o)*(x1(n)-x1o)
           swy = swy + uu(n)*(x2(n)-x2o)*(y1(n)-y1o)
           sw  = sw  + uu(n)*(x2(n)-x2o)
           szx = szx + uu(n)*(y2(n)-y2o)*(x1(n)-x1o)
           szy = szy + uu(n)*(y2(n)-y2o)*(y1(n)-y1o)
           sz  = sz  + uu(n)*(y2(n)-y2o)
           enddo
 
      dsxx = sx*sx - NNN*sxx
      dsyy = sy*sy - NNN*syy
      dsxy = sx*sy - NNN*sxy
 
      dlta = dsxx*dsyy - dsxy*dsxy
 
      if (dlta.eq.0) goto 999
 
      dswx = sw*sx - NNN*swx
      dswy = sw*sy - NNN*swy
      dszx = sz*sx - NNN*szx
      dszy = sz*sy - NNN*szy
 
      A =   (dswx*dsyy-dswy*dsxy)/dlta
      B =   (dswy*dsxx-dswx*dsxy)/dlta
      C =   (dszx*dsyy-dszy*dsxy)/dlta
      D =   (dszy*dsxx-dszx*dsxy)/dlta
 
 
      return
 
 999  continue
      A  = 1
      B  = 0
      C  = 0
      D  = 1
      x1o = 0
      y1o = 0
      x2o = 0
      y2o = 0
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/ROUTINES/TRANS/glob_fit6nrDPU_REJ.f"
c****                                         
c*********************************************
                                              
 
 
c---------------------------------------------------------------------
c
c    RMS2D:  report the fitting residual
c SIG_CLIP:  the sigma clipping threshold
c  RES_MIN:  never reject anything within this residual (where to stop)
c  RES_MAX:  never accept anything within this residual
c
c---------------------------------------------------------------------
 
      subroutine glob_fit6nrDPU_REJ(x1,y1,x2,y2,uu,Ns,Nu,           ! added Nu
     .                              A,B,C,D,x1o,y1o,x2o,y2o,
     .                              RMS2D,SIG_CLIP,RES_MIN,RES_MAX) ! added sev
      implicit none
 
      real*8 x1(*), y1(*)
      real*8 x2(*), y2(*)
      integer uu(*)
      integer Ns
      integer Nu
      real*8 A, B, C, D
      real*8 x1o, y1o, x2o, y2o
      real   RMS2D
      real   SIG_CLIP
      real   RES_MIN
      real   RES_MAX
 
      integer n, Us
      real dx, dy, dd
 
      real    dmax
      integer nmax
      real    drms, drmsX
 
 
      drms = 0.
 
    1 continue
      Us = 0
      do N = 1, Ns
         Us = Us + uu(N)
         enddo
      if (Us.le.3) goto 3
 
      call glob_fit6nrDPU(x1,y1,x2,y2,uu,Ns,
     .                    A,B,C,D,x1o,y1o,x2o,y2o)
 
      drms = 0.
      dmax = 0.
      nmax = 0
      do N = 1, Ns
         if (uu(N).eq.1) then
             dx = SNGL(x2(N)-x2o-A*(x1(N)-x1o)-B*(y1(N)-y1o))
             dy = SNGL(y2(N)-y2o-C*(x1(N)-x1o)-D*(y1(N)-y1o))
             dd = sqrt(dx**2+dy**2)
             if (dd.gt.dmax) then
                 dmax = dd
                 nmax = n
                 endif
             drms = drms + dd**2
             endif
         enddo
 
      drmsX  = drms - dmax**2
      drmsX  = sqrt(drmsX/(Us-2))
      drms   = sqrt(drms /(Us-1))
 
      if (dmax.le.RES_MIN) goto 3 ! everything is all within specs
 
      if (dmax.gt.RES_MAX) then ! one point is beyond specs
          uu(nmax) = 0
          goto 1
          endif
 
      if (dmax.gt.drmsX*SIG_CLIP) then ! sigma-clip this point
          uu(nmax) = 0
          goto 1
          endif
 
    3 continue
      call glob_fit6nrDPU(x1,y1,x2,y2,uu,Ns,
     .                    A,B,C,D,x1o,y1o,x2o,y2o)
 
      Nu = Us
      RMS2D = drms
 
 
      return
      end
 
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/PARLOR/islinux.f"
c****                                         
c*********************************************
                                              
 
c      program q
c      implicit none
c
c      logical islinux, islinux_var
c
c      islinux_var = islinux()
c
c      islinux_var = islinux()
c
c      stop
c      end
 
 
      logical function islinux()
      implicit none
 
      logical islinux_save
 
      integer i
      data i/0/
      data islinux_save/.true./
 
      common /islinux_/i,islinux_save
 
      byte    b(2)
      equivalence(i,b)
 
      if (i.eq.1) goto 1
 
      i = 1
      islinux_save = .false.
      if (b(1).eq.1) islinux_save = .true.
 
    1 islinux = islinux_save
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/notblank70.f"
c****                                         
c*********************************************
                                              
 
c----------------------------------------------
c
c
      logical function notblank70(string70)
      implicit none
 
      character*70 string70
 
      integer i
 
      notblank70 = .false.
      do i = 1, 70
         if (string70(i:i).ge.'a'.and.
     .       string70(i:i).le.'z') notblank70 = .true.
         if (string70(i:i).ge.'A'.and.
     .       string70(i:i).le.'Z') notblank70 = .true.
         if (string70(i:i).ge.'0'.and.
     .       string70(i:i).le.'9') notblank70 = .true.
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/writfits_b1.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine writfits_b1(FILE,pix,PXDIMX,PXDIMY)
      implicit none
 
      character*(*) FILE
      integer PXDIMX,PXDIMY
      byte pix(PXDIMX,PXDIMY)
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
 
      character*70 HDR(25)
      common/HDR/HDR
 
      character*80 FILEU
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
 
 
c     write(*,'(''ENTER   writfits_b1: '',80a)') FILEU
 
      open(10,file=FILEU,
     .     status='unknown',
     .     err   =900,
     .     recl  =2880,
     .     form  ='UNFORMATTED',
     .     access='DIRECT')
 
 
      i = 1
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  = T    '')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =  8  '')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)')
     .                          ' ''BYTE    '''
      write(buffc( 6*80+1: 7*80),'(''DATE    = '',11a)')
     .                          ' ''28/01/00'''
      write(buffc( 7*80+1: 8*80),'(''BSCALE  = '',i12)') 00001
      write(buffc( 8*80+1: 9*80),'(''BZERO   = '',i12)') 00000
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a20)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a20)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a20)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a20)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a20)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a20)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   = '',a20)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   = '',a20)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   = '',a20)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   = '',a20)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a20)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER = '',a20)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   = '',a20)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a20)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a20)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME = '',a20)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a20)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME= '',a20)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG = '',a20)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a20)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID= '',a20)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 = '',a20)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 = '',a20)') HDR(23)
      write(buffc(32*80+1:33*80),'(''VAFACTOR= '',a20)') HDR(24)
      write(buffc(33*80+1:34*80),'(''COMMENT   '')')
      write(buffc(34*80+1:35*80),'(''COMMENT   '')')
      write(buffc(35*80+1:36*80),'(''END       '')')
 
c     print*,buffc
      write(10,rec=i,iostat=ios) buffc
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  =   PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
 
c      print*,'nbyte1: ',nbyte1
c      print*,'nbyte2: ',nbyte1
c      print*,'    i1: ',i1
c      print*,'    i2: ',i2
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1) + 1
         np2 = (nbyteE-nbyte1)+ 1
        !if (i.lt.0010) print*,'i: ',i,i1,i2,np1,np2,nbyte0
        !if (i.gt.2900) print*,'i: ',i,i1,i2,np1,np2,nbyte0
         call pix2buff_b1(pix,buffb,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo
 
      close(10)
 
      return
 
 900  continue
      print*,'WRITFITS.f ERROR'
      stop
 
      end
 
      subroutine pix2buff_b1(pix,buffb,n1,nt)
      implicit none
      byte pix(1)
      byte buffb(2880)
 
      integer   n1,nt
 
      integer i
      integer npu
 
      byte ipval
 
      do i = 1, 2880
         npu = (1+(n1-1)*2)+i-1
         npu = (1+(n1-1)*1)+i-1
         if (npu.ge.1.and.npu.le.nt) then
            ipval = pix(npu)
            buffb(i) = ipval
            endif
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/writfits_r4.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------------
c
c this just writes a real*4 fits image
c
 
      subroutine writfits_r4(FILE,pix,PXDIMX,PXDIMY)
      implicit none
 
      character*(*) FILE
      integer PXDIMX,PXDIMY
      real    pix(PXDIMX,PXDIMY)
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffb,buffc)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer k, j
 
      character*200 FILEU
      character*70 HDR(25)
      common/HDR/HDR
 
      integer nonspace
      logical notblank70
 
      HDR(05) = "' pixel    '"
      HDR(06) = "' pixel    '"
 
      FILEU = FILE
      do i = 1, 195
         if (FILE(i:i+4).eq.'.fits') then
             FILEU = FILE(1:i+4)
             goto 1
             endif
         enddo
 
 
      print*,' '
      print*,'writfits_r4:    no .fits'
      print*,'      FILE:  ',FILE
      stop
    1 continue
 
      i = 1
 
      open(10,file=FILEU,status='unknown',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      do k = 00, 35
         write(buffc(k*80+1:k*80+80),'(80x)')
         enddo
 
      write(10,rec=i,iostat=ios) buffc
 
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  =                    T'')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =                  -32'')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',8x,i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',8x,i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',8x,i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)')
     .                          " 'REAL*4' "
      write(buffc(07*80+1:08*80),'(''COMMENT   '',a05)') '     '
      write(buffc(08*80+1:09*80),'(''COMMENT   '',a05)') '     '
      if (notblank70(HDR(01)))
     .write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a70)') HDR(01)
      if (notblank70(HDR(02)))
     .write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a70)') HDR(02)
      if (notblank70(HDR(03)))
     .write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a70)') HDR(03)
      if (notblank70(HDR(04)))
     .write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a70)') HDR(04)
      if (notblank70(HDR(05)))
     .write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a70)') HDR(05)
      if (notblank70(HDR(06)))
     .write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a70)') HDR(06)
      if (notblank70(HDR(07)))
     .write(buffc(15*80+1:16*80),'(''CD1_1   = '',a70)') HDR(07)
      if (notblank70(HDR(08)))
     .write(buffc(16*80+1:17*80),'(''CD1_2   = '',a70)') HDR(08)
      if (notblank70(HDR(09)))
     .write(buffc(17*80+1:18*80),'(''CD2_1   = '',a70)') HDR(09)
      if (notblank70(HDR(10)))
     .write(buffc(18*80+1:19*80),'(''CD2_2   = '',a70)') HDR(10)
      if (notblank70(HDR(11)))
     .write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a70)') HDR(11)
      if (notblank70(HDR(12)))
     .write(buffc(20*80+1:21*80),'(''PA_APER = '',a70)') HDR(12)
      if (notblank70(HDR(13)))
     .write(buffc(21*80+1:22*80),'(''PA_V3   = '',a70)') HDR(13)
      if (notblank70(HDR(14)))
     .write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a70)') HDR(14)
      if (notblank70(HDR(15)))
     .write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a70)') HDR(15)
      if (notblank70(HDR(16)))
     .write(buffc(24*80+1:25*80),'(''EXPTIME = '',a70)') HDR(16)
      if (notblank70(HDR(17)))
     .write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a70)') HDR(17)
      if (notblank70(HDR(18)))
     .write(buffc(26*80+1:27*80),'(''TARGNAME= '',a70)') HDR(18)
      if (notblank70(HDR(19)))
     .write(buffc(27*80+1:28*80),'(''RA_TARG = '',a70)') HDR(19)
      if (notblank70(HDR(20)))
     .write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a70)') HDR(20)
      if (notblank70(HDR(21)))
     .write(buffc(29*80+1:30*80),'(''PROPOSID= '',a70)') HDR(21)
      if (notblank70(HDR(22)))
     .write(buffc(30*80+1:31*80),'(''FILTER1 = '',a70)') HDR(22)
      if (notblank70(HDR(23)))
     .write(buffc(31*80+1:32*80),'(''FILTER2 = '',a70)') HDR(23)
      if (notblank70(HDR(24)))
     .write(buffc(33*80+1:34*80),'(''VAFACTOR= '',a70)') HDR(24)
      if (notblank70(HDR(25)))
     .write(buffc(32*80+1:33*80),'(''CCDGAIN = '',a70)') HDR(25)
      write(buffc(33*80+1:34*80),'(''COMMENT   '',a05)') '     '
      write(buffc(34*80+1:35*80),'(''COMMENT   '',a05)') '     '
      write(buffc(35*80+1:36*80),'(''END       '')')
 
      do k = 00, 34
         nonspace = 0
         do j = 10, 80
            if (buffc(k*80+j:k*80+j).ne.' '.and.
     .          iachar(buffc(k*80+j:k*80+j)).ne.0) nonspace = nonspace+1
            enddo
         if (nonspace.eq.0) then
             do j = 1, 80
                buffc(k*80+j:k*80+j) = ' '
                enddo
             endif
         enddo
 
      write(10,rec=i,iostat=ios) buffc
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 4*PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         call pix2buff_r4(buffb,pix,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo
 
      close(10)
 
      return
 
 900  continue
      print*,'writfits_r4.f ERROR'
      print*,'   FILEU: ',FILEU
      stop
 
      end
 
 
 
c-------------------------------------------
c
c
c
      logical function notblank70_(string70)
      implicit none
 
      character*70 string70
 
      integer i
 
      notblank70_ = .false.
      do i = 1, 70
         if (string70(i:i).ge.'a'.and.
     .       string70(i:i).le.'z') notblank70_ = .true.
         if (string70(i:i).ge.'A'.and.
     .       string70(i:i).le.'Z') notblank70_ = .true.
         if (string70(i:i).ge.'0'.and.
     .       string70(i:i).le.'9') notblank70_ = .true.
         enddo
 
      return
      end
 
 
 
c-------------------------------------------------------
c
c
      subroutine pix2buff_r4(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      real*4 pix(*)
      integer n1,nt
 
      byte b(4)
      real*4 r
      equivalence(r,b)
 
      integer i, npu, nbu
      logical islinux
 
      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) r = pix(npu)
         if (.not.islinux()) then
            buff(nbu+1) = b(1)
            buff(nbu+2) = b(2)
            buff(nbu+3) = b(3)
            buff(nbu+4) = b(4)
            endif
         if (islinux()) then
            buff(nbu+1) = b(4)
            buff(nbu+2) = b(3)
            buff(nbu+3) = b(2)
            buff(nbu+4) = b(1)
            endif
         enddo
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/query_hdr.f"
c****                                         
c*********************************************
                                              
 
c---------------------------------------------------
c
c This routine checks to see if the header of an image
c contains a particular keyword. If it does, it returns
c its value-string in STREAMX.
c
 
      subroutine query_hdr(filename,FIELDX,streamx)
      implicit none
 
      character*200 filename
 
      character*8  field
      character*20 stream
 
      character*8  fieldx
      character*20 streamx
 
      integer  i
      integer ios, k
 
      character*2880 buff
 
      integer nread
      integer NEXTEND
 
c-----------------------------------------------
 
      !print*,'query_hdr...',FILENAME
      !print*,'    fieldx: ',fieldx
      streamx = ' '
 
      close(10)
      open(10,file=FILENAME,status='old',iostat=ios,
     .  err=900,recl=2880,form='UNFORMATTED',
     .  access='DIRECT')
 
      streamx = 'NULL'
      i = 0
      NREAD   = 0
      NEXTEND = 0
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios,err=901) buff
      if (ios.lt.0) goto 900
      do k = 0, 35, 1
         field  = buff(k*80+01:k*80+08)
         stream = buff(k*80+11:k*80+31)
         !print*,i,k,field,stream
         if (field.eq.'NEXTEND ') read(stream,*) NEXTEND
         if (field.eq.fieldx) streamx = stream(1:20)
         if (field.eq.'END     ') goto 101
         continue
         enddo
      goto 100
 101  continue
      NREAD = NREAD + 1
      if (NREAD.le.1.and.NEXTEND.gt.1) goto 100
      close(10)
      !print*,'    streamx: ',streamx
 
      return
 
  900 continue
      print*,' '
      print*,'query_hdr() FILE OPEN ERROR EXIT. '
      print*,' '
      print*,'ONE OF THE IMAGES WAS NOT IN STANDARD'
      print*,'HST FITS FORMAT.'
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      write(*,'(''      FIELDX: '',a08)') FIELDX
      write(*,'(''      IOSTAT: '',i4 )') IOS
      print*,' '
      stop
 
  901 continue
      print*,' '
      print*,'query_hdr() FILE READ ERROR EXIT. '
      print*,' '
      print*,'ONE OF THE IMAGES WAS NOT IN STANDARD'
      print*,'HST FITS FORMAT.'
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      write(*,'(''      FIELDX: '',a08)') FIELDX
      write(*,'(''      IOSTAT: '',i4 )') IOS
      print*,' '
      stop
 
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/readfits_r4e.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------------------------
c
c reads an r4 fits image, with extensions (reads in one extension)
c
 
      subroutine readfits_r4e(FILE,pix,NDIMX,NDIMY,NEXTENU)
      implicit none
 
      character*(*) FILE
      integer NDIMX,NDIMY
      real    pix(NDIMX,NDIMY)
      integer NEXTENU
 
      character*200 FILEU
      character*70 INFO(10)
      common / fitsinfo / INFO
 
      integer naxes
      integer laxis(3)
      common/laxis3_/laxis
 
      integer NXF
      integer NYF
 
      character*8  field
      character*40 stream
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
c     integer j
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
      character*70 HDR(25)
      common/HDR/HDR
 
      logical DIAG
      data DIAG /.false./
      integer NEND
      integer iend
 
      if (DIAG) then
         print*,'enter readfits_r4e...... '
         print*,'FILEi: ',FILE(1:60)
         endif
 
      FILEU = FILE
      iend = 0
      do i = 196,1,-1
         if (DIAG) print*,i,iend,FILE(i:i+4)
         if (FILE(i:i+4).eq.'.fits') iend = i+4
         enddo
      if (DIAG) then
         print*,'iend: ',iend
         endif
      if (iend.eq.0) then
          print*,'NO .fits in FILENAME'
          print*,'FILE: ',FILE
          stop
          endif
      FILEU = FILE(1:iend)
 
      if (DIAG) then
         print*,'purge...'
         print*,'FILEu: ',FILEU
         endif
 
      open(10,file=FILEU,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      if (DIAG) print*,'...opened ',FILEU
 
      naxes = -1
      laxis(1) = 1
      laxis(2) = 1
      laxis(3) = 1
      nextend = 0
 
      do i = 1, 10
         INFO(i) = ' '
         enddo
 
      BSCALE = 1.0
      BZERO  = 0.0
 
      NEND = 0
      i = 0
      nread = 0
 100  continue
      i = i + 1
      if (DIAG) print*,'READREC: ',i
      read(10,rec=i,iostat=ios) buffc
      do k = 0, 35, 1
         if (DIAG) write(*,'(i4,1x,i4,1x,a80)')
     .                   i,k,buffc(k*80+1:k*80+80)
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+11:k*80+51)
         if (field.eq.'NAXIS   ') read(stream,*) naxes
         if (field.eq.'NAXIS1  ') read(stream,*) laxis(1)
         if (field.eq.'NAXIS2  ') read(stream,*) laxis(2)
         if (field.eq.'NAXIS3  ') read(stream,*) laxis(3)
         if (field.eq.'NEXTEND ') read(stream,*) nextend
         if (field.eq.'BITPIX  ') read(stream,*) bitpix
         if (field.eq.'BSCALE  ') read(stream,*) bscale
         if (field.eq.'BZERO   ') read(stream,*) bzero
 
         if (field.eq.'EXPTIME ') INFO(1) = stream
         if (field.eq.'FILTER  ') INFO(2) = stream
         if (field.eq.'FILTNAM1') INFO(2) = stream
         if (field.eq.'FILENAME') INFO(3) = stream
         if (field.eq.'DATE-OBS') INFO(4) = stream
         if (field.eq.'TIME-OBS') INFO(5) = stream
         if (field.eq.'DEC_TARG') INFO(6) = stream
         if (field.eq.'RA_TARG ') INFO(7) = stream
         if (field.eq.'PA_V3   ') INFO(8) = stream
         if (field.eq.'PROPOSID') INFO(9) = stream
 
         if (field.eq.'CRPIX1  ') HDR(01) = stream
         if (field.eq.'CRPIX2  ') HDR(02) = stream
         if (field.eq.'CRVAL1  ') HDR(03) = stream
         if (field.eq.'CRVAL2  ') HDR(04) = stream
         if (field.eq.'CTYPE1  ') HDR(05) = stream
         if (field.eq.'CTYPE2  ') HDR(06) = stream
         if (field.eq.'CD1_1   ') HDR(07) = stream
         if (field.eq.'CD1_2   ') HDR(08) = stream
         if (field.eq.'CD2_1   ') HDR(09) = stream
         if (field.eq.'CD2_2   ') HDR(10) = stream
         if (field.eq.'ORIENTAT') HDR(11) = stream
         if (field.eq.'PA_APER ') HDR(12) = stream
         if (field.eq.'PA_V3   ') HDR(13) = stream
         if (field.eq.'DATE-OBS') HDR(14) = stream
         if (field.eq.'TIME-OBS') HDR(15) = stream
         if (field.eq.'EXPTIME ') HDR(16) = stream
         if (field.eq.'ROOTNAME') HDR(17) = stream
         if (field.eq.'TARGNAME') HDR(18) = stream
         if (field.eq.'RA_TARG ') HDR(19) = stream
         if (field.eq.'DEC_TARG') HDR(20) = stream
         if (field.eq.'PROPOSID') HDR(21) = stream
         if (field.eq.'FILTER ')  HDR(22) = stream
         if (field.eq.'FILTER1 ') HDR(22) = stream
         if (field.eq.'FILTER2 ') HDR(23) = stream
         if (field.eq.'VAFACTOR') HDR(24) = stream
 
         if (field.eq.'END     ') then
           !print*,' ---> NEND: ',NEND,NEXTENU
            if (NEND.eq.NEXTENU) goto 101
            NEND = NEND + 1
            endif
         enddo
      goto 100
 101  continue
 
      nread = nread + 1
      if (DIAG) then
         print*,' '
         print*,'----------------------------------------'
         print*,'  NREAD: ',nread
         print*,'NEXTEND: ',nextend
         print*,'  NAXIS: ',naxes
         print*,'  LAXIS: ',laxis(1),laxis(2),laxis(3)
         print*,' BITPIX: ',bitpix
         print*,' BSCALE: ',bscale
         print*,'  BZERO: ',bzero
         print*,'  NDIMX: ',NDIMX
         print*,'  NDIMY: ',NDIMY
         print*,' '
         endif
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      NXF = laxis(1)
      NYF = laxis(2)
 
      nbper  = 4*laxis(1)*laxis(2)
      npt    = laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      if (BITPIX.eq.-64) then
         print*,'from readfits_r4e... '
         print*,'     BITPIX: ',BITPIX
         print*,'call readfits_r8e... '
         print*,'                   '
         call readfits_r8e_r4(FILE,pix,NDIMX,NDIMY,NEXTENU)
         return
         endif
 
      if (BITPIX.ne.-32) then
         print*,'readfits_r4e...'
         print*,'BITPIX: ',BITPIX
         print*,'prob'
         stop
         endif
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
c        call buff2pix_r4(buffb,pix,np1,npt)
         call buff2pix_r4_edge(buffb,pix,np1,npt,
     .                         NDIMX,NDIMY,laxis(1),laxis(2))
         if (DIAG) write(*,1115) i,np1,np2,npt
 1115    format(1x,i8,1x,i10,1x,i10,1x,i10)
         enddo
 
      close(10)
      if (DIAG) print*,'...closed ',FILEU
      return
 
 900  continue
      print*,' '
      print*,'readfits_r4e ERROR'
      write(*,'(''NEXTU: '',i4  )') NEXTENU
      write(*,'(''FILEI: '',a200)') FILE
      write(*,'('' IEND: '',i4  )') IEND
      write(*,'(''FILEU: '',a200)') FILEU
      print*,' '
      stop
 
      end
 
 
c------------------------------------------------------
c
c
      subroutine buff2pix_r4_edge(buff,pix,n1,nt,
     .                            NXP,NYP,NXF,NYF)
      implicit none
c     character buff(2880)
      byte buff(2880)
      integer NXP,NYP
      real pix(NXP,NYP)
      integer n1,nt
      integer NXF,NYF
 
      real pbuff(720)
      common /sneaky/pbuff
      integer i
      integer npu, nbu
      integer NX, NY
 
      byte b(4)
      real r
      equivalence(r,b)
      logical islinux
 
      if (NYF.eq.0) continue ! just to "use" the argument
 
      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) then
            NX = npu - (npu-1)/NXF*NXF
            NY = 1 + (npu-1)/NXF
            if (NX.le.NXP.and.NY.le.NYP) then
               if (.not.islinux()) then
                  b(1) = buff(nbu+1)
                  b(2) = buff(nbu+2)
                  b(3) = buff(nbu+3)
                  b(4) = buff(nbu+4)
                  endif
               if (islinux()) then
                  b(4) = buff(nbu+1)
                  b(3) = buff(nbu+2)
                  b(2) = buff(nbu+3)
                  b(1) = buff(nbu+4)
                  endif
               if (npu.ge.1.and.npu.le.nt) pix(NX,NY) = r
               endif
            endif
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/readfits_r8e_r4.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------------------------
c
c reads an r4 fits image, with extensions (reads in one extension)
c
 
      subroutine readfits_r8e_r4(FILE,pix,NDIMX,NDIMY,NEXTENU)
      implicit none
 
      character*(*) FILE
      integer NDIMX,NDIMY
      real*4  pix(NDIMX,NDIMY)
      integer NEXTENU
 
      character*199 FILEU
      character*70 INFO(10)
      common / fitsinfo / INFO
 
      integer naxes
      integer laxis(3)
      common/laxis3_/laxis
 
      integer NXF
      integer NYF
 
      character*8  field
      character*40 stream
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
      integer j
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
      character*70 HDR(25)
      common/HDR/HDR
 
      logical DIAG
      data DIAG /.false./
      integer NEND
      integer iend
 
      if (DIAG) then
         print*,'enter readfits_r4e...... '
         print*,'FILEi: ',FILE(1:60)
         endif
 
      FILEU = FILE
      iend = 0
      do i = 196,1,-1
         if (DIAG) print*,i,iend,FILE(i:i+4)
         if (FILE(i:i+4).eq.'.fits') iend = i+4
         enddo
      if (DIAG) then
         print*,'iend: ',iend
         endif
      if (iend.eq.0) stop 'NO .fits in FILENAME'
      FILEU = FILE(1:iend)
 
      if (DIAG) then
         print*,'purge...'
         print*,'FILEu: ',FILEU
         endif
 
      open(10,file=FILEU,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      if (DIAG) print*,'...opened ',FILEU
 
      naxes = -1
      laxis(1) = 1
      laxis(2) = 1
      laxis(3) = 1
      nextend = 0
 
      do i = 1, 10
         INFO(i) = ' '
         enddo
 
      BSCALE = 1.0
      BZERO  = 0.0
 
      NEND = 0
      i = 0
      nread = 0
 100  continue
      i = i + 1
      if (DIAG) print*,'READREC: ',i
      read(10,rec=i,iostat=ios) buffc
      do k = 0, 35, 1
         if (DIAG) write(*,'(i4,1x,i4,1x,a80)')
     .                   i,k,buffc(k*80+1:k*80+80)
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+11:k*80+51)
         if (field.eq.'NAXIS   ') read(stream,*) naxes
         if (field.eq.'NAXIS1  ') read(stream,*) laxis(1)
         if (field.eq.'NAXIS2  ') read(stream,*) laxis(2)
         if (field.eq.'NAXIS3  ') read(stream,*) laxis(3)
         if (field.eq.'NEXTEND ') read(stream,*) nextend
         if (field.eq.'BITPIX  ') read(stream,*) bitpix
         if (field.eq.'BSCALE  ') read(stream,*) bscale
         if (field.eq.'BZERO   ') read(stream,*) bzero
 
         if (field.eq.'EXPTIME ') INFO(1) = stream
         if (field.eq.'FILTER  ') INFO(2) = stream
         if (field.eq.'FILTNAM1') INFO(2) = stream
         if (field.eq.'FILENAME') INFO(3) = stream
         if (field.eq.'DATE-OBS') INFO(4) = stream
         if (field.eq.'TIME-OBS') INFO(5) = stream
         if (field.eq.'DEC_TARG') INFO(6) = stream
         if (field.eq.'RA_TARG ') INFO(7) = stream
         if (field.eq.'PA_V3   ') INFO(8) = stream
         if (field.eq.'PROPOSID') INFO(9) = stream
 
         if (field.eq.'CRPIX1  ') HDR(01) = stream
         if (field.eq.'CRPIX2  ') HDR(02) = stream
         if (field.eq.'CRVAL1  ') HDR(03) = stream
         if (field.eq.'CRVAL2  ') HDR(04) = stream
         if (field.eq.'CTYPE1  ') HDR(05) = stream
         if (field.eq.'CTYPE2  ') HDR(06) = stream
         if (field.eq.'CD1_1   ') HDR(07) = stream
         if (field.eq.'CD1_2   ') HDR(08) = stream
         if (field.eq.'CD2_1   ') HDR(09) = stream
         if (field.eq.'CD2_2   ') HDR(10) = stream
         if (field.eq.'ORIENTAT') HDR(11) = stream
         if (field.eq.'PA_APER ') HDR(12) = stream
         if (field.eq.'PA_V3   ') HDR(13) = stream
         if (field.eq.'DATE-OBS') HDR(14) = stream
         if (field.eq.'TIME-OBS') HDR(15) = stream
         if (field.eq.'EXPTIME ') HDR(16) = stream
         if (field.eq.'ROOTNAME') HDR(17) = stream
         if (field.eq.'TARGNAME') HDR(18) = stream
         if (field.eq.'RA_TARG ') HDR(19) = stream
         if (field.eq.'DEC_TARG') HDR(20) = stream
         if (field.eq.'PROPOSID') HDR(21) = stream
         if (field.eq.'FILTER ')  HDR(22) = stream
         if (field.eq.'FILTER1 ') HDR(22) = stream
         if (field.eq.'FILTER2 ') HDR(23) = stream
         if (field.eq.'VAFACTOR') HDR(24) = stream
 
         if (field.eq.'END     ') then
           !print*,' ---> NEND: ',NEND,NEXTENU
            if (NEND.eq.NEXTENU) goto 101
            NEND = NEND + 1
            endif
         enddo
      goto 100
 101  continue
 
      nread = nread + 1
      if (DIAG) then
         print*,' '
         print*,'----------------------------------------'
         print*,'  NREAD: ',nread
         print*,'NEXTEND: ',nextend
         print*,'  NAXIS: ',naxes
         print*,'  LAXIS: ',laxis(1),laxis(2),laxis(3)
         print*,' BITPIX: ',bitpix
         print*,' BSCALE: ',bscale
         print*,'  BZERO: ',bzero
         print*,'  NDIMX: ',NDIMX
         print*,'  NDIMY: ',NDIMY
         print*,' '
         endif
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      NXF = laxis(1)
      NYF = laxis(2)
 
      nbper  = 8*laxis(1)*laxis(2)
      npt    = laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      if (BITPIX.ne.-64) then
         print*,'readfits_r8e...'
         print*,'BITPIX: ',BITPIX
         print*,'prob'
         stop
         endif
 
      if (DIAG) then
         print*,'          '
         print*,'do i ---  ',i1,i2
         print*,'   nbper: ',nbper
         print*,'     npt: ',npt
         print*,'          '
         endif
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/8 + 1
         np2 = (nbyteE-nbyte1)/8 + 1
         call buff2pix_r8_edge(buffb,pix,np1,npt,
     .                         NDIMX,NDIMY,laxis(1),laxis(2))
         if (DIAG) write(*,1115) i,np1,np2,npt
 1115    format(1x,i8,1x,i10,1x,i10,1x,i10)
         enddo
 
      if (DIAG) then
          print*,' '
          do j = 10, 01, -1
             write(*,'(10x,i2,1x,10f8.2)') j, (pix(i,j),i=1,10)
             enddo
          print*,' '
          endif
 
      close(10)
      if (DIAG) print*,'...closed ',FILEU
      return
 
 900  continue
      print*,'READFITS ERROR'
      stop
 
      end
 
 
c------------------------------------------------------
c
c
      subroutine buff2pix_r8_edge(buff,pix,n1,nt,
     .                            NXP,NYP,NXF,NYF)
      implicit none
      byte buff(2880)
      integer NXP,NYP
      real*4  pix(NXP,NYP)
      integer n1,nt
      integer NXF,NYF
 
c     real pbuff(360)
      integer i
      integer npu, nbu
      integer NX, NY
 
      byte b(8)
      real*8 r8
      equivalence(r8,b)
      logical islinux
 
      if (NYF.eq.0) continue
 
      do i = 1, 360
         npu = n1+i-1
         nbu = (i-1)*8
         if (npu.ge.1.and.npu.le.nt) then
            NX = npu - (npu-1)/NXF*NXF
            NY = 1 + (npu-1)/NXF
            if (NX.le.NXP.and.NY.le.NYP) then
               if (.not.islinux()) then
                  b(1) = buff(nbu+1)
                  b(2) = buff(nbu+2)
                  b(3) = buff(nbu+3)
                  b(4) = buff(nbu+4)
                  b(5) = buff(nbu+5)
                  b(6) = buff(nbu+6)
                  b(7) = buff(nbu+7)
                  b(8) = buff(nbu+8)
                  endif
               if (islinux()) then
                  b(8) = buff(nbu+1)
                  b(7) = buff(nbu+2)
                  b(6) = buff(nbu+3)
                  b(5) = buff(nbu+4)
                  b(4) = buff(nbu+5)
                  b(3) = buff(nbu+6)
                  b(2) = buff(nbu+7)
                  b(1) = buff(nbu+8)
                  endif
               if (npu.ge.1.and.npu.le.nt) pix(NX,NY) = SNGL(r8)
               endif
            endif
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/writfits_r4_WCSLTV.f"
c****                                         
c*********************************************
                                              
 
 
c-----------------------------------------------------
c
c this just writes a real*4 fits image
c
 
      subroutine writfits_r4_WCSLTV(FILE,pix,PXDIMX,PXDIMY,
     .                              CRPIX1,CRPIX2,CRVAL1,CRVAL2,
     .                              CD1_1, CD1_2, CD2_1, CD2_2,
     .                              LTV1, LTV2,
     .                              LTM1_1,LTM1_2,LTM2_1,LTM2_2)
      implicit none
 
      character*(*) FILE
      integer PXDIMX,PXDIMY
      real    pix(PXDIMX,PXDIMY)
      real*8  CRPIX1,CRPIX2,CRVAL1,CRVAL2
      real*8  CD1_1, CD1_2, CD2_1, CD2_2
      real*8  LTV1,  LTV2
      real*8  LTM1_1,LTM1_2,LTM2_1,LTM2_2
 
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffb,buffc)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer k
 
      character*80 FILEU
      character*70 HDR(25)
      common/HDR/HDR
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
      open(10,file=FILEU,status='unknown',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      do k = 00, 34
         write(buffc(k*80+1:k*80+80),'(80x)')
         enddo
 
      write(10,rec=i,iostat=ios) buffc
 
      write(buffc(00*80+1:01*80),'(''SIMPLE  =                    T'')')
      write(buffc(01*80+1:02*80),'(''BITPIX  =                  -32'')')
      write(buffc(02*80+1:03*80),'(''NAXIS   = '',8x,i12)') 2
      write(buffc(03*80+1:04*80),'(''NAXIS1  = '',8x,i12)') PXDIMX
      write(buffc(04*80+1:05*80),'(''NAXIS2  = '',8x,i12)') PXDIMY
      write(buffc(05*80+1:06*80),'(''DATATYPE= '',9a)')  " 'REAL*4' "
      write(buffc(06*80+1:07*80),'(''COMMENT   '',a05)') '     '
      write(buffc(07*80+1:08*80),'(''CTYPE1  = '',a10)') "'RA---TAN'"
      write(buffc(08*80+1:09*80),'(''CTYPE2  = '',a10)') "'DEC--TAN'"
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',f20.10)') CRPIX1
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',f20.10)') CRPIX2
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',f20.10)') CRVAL1
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',f20.10)') CRVAL2
      write(buffc(13*80+1:14*80),'(''CD1_1   = '',f20.10)') CD1_1
      write(buffc(14*80+1:15*80),'(''CD1_2   = '',f20.10)') CD1_2
      write(buffc(15*80+1:16*80),'(''CD2_1   = '',f20.10)') CD2_1
      write(buffc(16*80+1:17*80),'(''CD2_2   = '',f20.10)') CD2_2
      write(buffc(17*80+1:18*80),'(''COMMENT   '',a05)') '     '
      write(buffc(18*80+1:19*80),'(''LTV1    = '',f20.10)') LTV1
      write(buffc(19*80+1:20*80),'(''LTV2    = '',f20.10)') LTV2
      write(buffc(20*80+1:21*80),'(''LTM1_1  = '',f20.10)') LTM1_1
      write(buffc(21*80+1:22*80),'(''LTM1_2  = '',f20.10)') LTM1_2
      write(buffc(22*80+1:23*80),'(''LTM2_1  = '',f20.10)') LTM2_1
      write(buffc(23*80+1:24*80),'(''LTM2_2  = '',f20.10)') LTM2_2
      write(buffc(24*80+1:25*80),'(''COMMENT   '',a05)') '     '
      write(buffc(25*80+1:26*80),'(''END       '')')
 
      write(10,rec=i,iostat=ios) buffc
 
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 4*PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         call pix2buff_r4_WCS(buffb,pix,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo
 
      close(10)
 
      return
 
 900  continue
      print*,'writfits_r4.f ERROR'
      print*,'   FILEU: ',FILEU
      stop
 
      end
 
c-------------------------------------------------------
c
c
      subroutine pix2buff_r4_WCS(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      real*4 pix(*)
      integer n1,nt
 
      byte b(4)
      real*4 r
      equivalence(r,b)
 
      integer i, npu, nbu
 
      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) r = pix(npu)
         if (.not.(_LINUX_)) then
            buff(nbu+1) = b(1)
            buff(nbu+2) = b(2)
            buff(nbu+3) = b(3)
            buff(nbu+4) = b(4)
            endif
         if ((_LINUX_)) then
            buff(nbu+1) = b(4)
            buff(nbu+2) = b(3)
            buff(nbu+3) = b(2)
            buff(nbu+4) = b(1)
            endif
         enddo
 
      return
      end
 
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/extract_stdgdc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
c
c
      subroutine extract_stdgdc(xr,yr,xc,yc,mc,FILE_STDGC,SENSE)
      implicit none
 
      real*8 xr, yr
      real*8 xc, yc, mc
      character*(*) FILE_STDGC
      integer SENSE
 
      save
 
      logical FIRST
      data    FIRST/.true./
 
      character*200 FILE_STDGCu
      character*200 FILE_STDGCc
 
      integer   NDIM_XGC, NDIM_YGC
      integer*4 XGC_0, NDIM_XCG
      integer*4 YGC_0, NDIM_YCG
 
      integer i, ios
      integer k
      character*08   FIELD
      character*20   STREAM
      character*2880 buffc
 
      real*8  fx,   fy
      real*8  xcu,  ycu
      integer ipix, jpix
 
      integer*4, dimension(:,:), allocatable :: xgc
      integer*4, dimension(:,:), allocatable :: ygc
      integer*2, dimension(:,:), allocatable :: mgc
      integer*4, dimension(:,:), allocatable :: xcg
      integer*4, dimension(:,:), allocatable :: ycg
 
      if (FILE_STDGC(1:4).eq.'NONE'.or.
     .    FILE_STDGC(1:4).eq.'NULL') then
         if (SENSE.eq.0) then
            mc = 0.
            return
            endif
         if (SENSE.eq.1) then ! forward...
            xc = xr
            yc = yr
            return
            endif
         if (SENSE.eq.-1) then ! forward...
            xr = xc
            yr = yc
            return
            endif
         endif
 
      if (FIRST) FILE_STDGCc = ' '
 
      FILE_STDGCu = FILE_STDGC
      do k = 196,2,-1
         if (FILE_STDGCu(k:k+4).eq.'.fits')
     .       FILE_STDGCu = FILE_STDGCu(1:k+4)
         enddo
 
      if (FILE_STDGCc.ne.FILE_STDGCu) then
          FILE_STDGCc = FILE_STDGCu
c         print*,'EXTRACT_STDGC --- OPEN(10): ',FILE_STDGCc
          open(10,file=FILE_STDGCc,
     .         status='old',
     .         err   =900,
     .         recl  =2880,
     .         form  ='UNFORMATTED',
     .         access='DIRECT')
          i = 1
          read(10,rec=i,iostat=ios) buffc
          do k = 0, 35
             FIELD  = buffc(k*80+01:k*80+08)
             STREAM = buffc(k*80+11:k*80+30)
             if (FIELD.eq.'NDIM_XGC') read(STREAM,*) NDIM_XGC
             if (FIELD.eq.'NDIM_YGC') read(STREAM,*) NDIM_YGC
             if (FIELD.eq.'NDIM_XCG') read(STREAM,*) NDIM_XCG
             if (FIELD.eq.'NDIM_YCG') read(STREAM,*) NDIM_YCG
             if (FIELD.eq.'XGC_0   ') read(STREAM,*) XGC_0
             if (FIELD.eq.'YGC_0   ') read(STREAM,*) YGC_0
             enddo
          if (.not.FIRST) deallocate(xgc,ygc,mgc,xcg,ycg)
          allocate(xgc(NDIM_XGC,NDIM_YGC))
          allocate(ygc(NDIM_XGC,NDIM_YGC))
          allocate(mgc(NDIM_XGC,NDIM_YGC))
          allocate(xcg(NDIM_XCG,NDIM_YCG))
          allocate(ycg(NDIM_XCG,NDIM_YCG))
          i = 2
          call getext_i4(xgc,NDIM_XGC,NDIM_YGC,i)
          call getext_i4(ygc,NDIM_XGC,NDIM_YGC,i)
          call getext_i2(mgc,NDIM_XGC,NDIM_YGC,i)
          call getext_i4(xcg,NDIM_XCG,NDIM_YCG,i)
          call getext_i4(ycg,NDIM_XCG,NDIM_YCG,i)
          close(10)
          FIRST = .false.
          endif
 
      if (SENSE.eq.1) then ! forward...
          ipix = int(xr)
          jpix = int(yr)
          xc = -9999
          yc = -9999
 
          if (ipix.lt.       1-1) return
          if (ipix.gt.NDIM_XGC+1) return
          if (jpix.lt.       1-1) return
          if (jpix.gt.NDIM_YGC+1) return
 
          if (ipix.lt.1) ipix = 1
          if (jpix.lt.1) jpix = 1
          if (ipix.gt.NDIM_XGC-1) ipix = NDIM_XGC-1
          if (jpix.gt.NDIM_YGC-1) jpix = NDIM_YGC-1
 
          fx = xr-ipix
          fy = yr-jpix
          xc = (1-fx)*(1-fy)*xgc(ipix  ,jpix  )/1d5
     .       + (1-fx)*( fy )*xgc(ipix  ,jpix+1)/1d5
     .       + ( fx )*(1-fy)*xgc(ipix+1,jpix  )/1d5
     .       + ( fx )*( fy )*xgc(ipix+1,jpix+1)/1d5
          yc = (1-fx)*(1-fy)*ygc(ipix  ,jpix  )/1d5
     .       + (1-fx)*( fy )*ygc(ipix  ,jpix+1)/1d5
     .       + ( fx )*(1-fy)*ygc(ipix+1,jpix  )/1d5
     .       + ( fx )*( fy )*ygc(ipix+1,jpix+1)/1d5
          return
          endif
 
      if (SENSE.eq.0) then ! find the zeropoint...
          ipix = int(xr)
          jpix = int(yr)
          mc = 0.000
          if (ipix.lt.         1) return
          if (ipix.gt.NDIM_XGC-1) return
          if (jpix.lt.         1) return
          if (jpix.gt.NDIM_YGC-1) return
          fx = xr-ipix
          fy = yr-jpix
          mc = (1-fx)*(1-fy)*mgc(ipix  ,jpix  )/1d4
     .       + (1-fx)*( fy )*mgc(ipix  ,jpix+1)/1d4
     .       + ( fx )*(1-fy)*mgc(ipix+1,jpix  )/1d4
     .       + ( fx )*( fy )*mgc(ipix+1,jpix+1)/1d4
          return
          endif
 
      if (SENSE.eq.-1) then ! reverse
          xcu =  1 + (xc-XGC_0)
          ycu =  1 + (yc-YGC_0)
          ipix = int(xcu)
          jpix = int(ycu)
          xr = -1
          yr = -1
          if (ipix.lt.         1) return
          if (ipix.gt.NDIM_XCG-1) return
          if (jpix.lt.         1) return
          if (jpix.gt.NDIM_YCG-1) return
          if (xcg(ipix  ,jpix  ).le.-2000000000) return
          if (xcg(ipix+1,jpix  ).le.-2000000000) return
          if (xcg(ipix  ,jpix+1).le.-2000000000) return
          if (xcg(ipix+1,jpix+1).le.-2000000000) return
           fx = xcu-ipix
           fy = ycu-jpix
           xr = (1-fx)*(1-fy)*xcg(ipix  ,jpix  )/1d5
     .        + (1-fx)*( fy )*xcg(ipix  ,jpix+1)/1d5
     .        + ( fx )*(1-fy)*xcg(ipix+1,jpix  )/1d5
     .        + ( fx )*( fy )*xcg(ipix+1,jpix+1)/1d5
           yr = (1-fx)*(1-fy)*ycg(ipix  ,jpix  )/1d5
     .        + (1-fx)*( fy )*ycg(ipix  ,jpix+1)/1d5
     .        + ( fx )*(1-fy)*ycg(ipix+1,jpix  )/1d5
     .        + ( fx )*( fy )*ycg(ipix+1,jpix+1)/1d5
           return
           endif
 
      print*,'ILLEGAL SENSE',SENSE
      stop
 
  900 continue
      print*,'file open error'
      print*,'  FILE_STDGCc: ',FILE_STDGCc
      print*,'  FILE_SGDGC : ',FILE_STDGC
      stop
      end
 
 
c--------------------------------------------------------------------
c
c
      subroutine getext_i4(pix,NDIMX,NDIMY,i)
      implicit none
 
      integer NDIMX, NDIMY
      integer*4 PIX(NDIMX,NDIMY)
      integer i
 
      integer o
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer ifirst, i1, i2, ios
      integer np1, np2, npt
      integer nbu, npu
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      byte      b(4)
      integer*4 ii
      equivalence(ii,b)
 
      integer ipix, jpix
 
      logical islinux
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 4*NDIMX*NDIMY
      npt    =   NDIMX*NDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         read(10,rec=i,iostat=ios) buffc
         do o = 1, 720
            npu = np1+o-1
            nbu = (o-1)*4
            ii  = 0
            ipix = npu - (npu-1)/NDIMX*NDIMX
            jpix =   1 + (npu-1)/NDIMX
            if (ipix.ge.1.and.ipix.le.NDIMX.and.
     .          jpix.ge.1.and.jpix.le.NDIMY) then
                if (.not.(islinux())) then
                   b(1) = buffb(nbu+1)
                   b(2) = buffb(nbu+2)
                   b(3) = buffb(nbu+3)
                   b(4) = buffb(nbu+4)
                   endif
                if ((islinux())) then
                   b(4) = buffb(nbu+1)
                   b(3) = buffb(nbu+2)
                   b(2) = buffb(nbu+3)
                   b(1) = buffb(nbu+4)
                   endif
                pix(ipix,jpix) = ii
                endif
            enddo!o
         enddo!i
 
      i = i2 + 1
 
      return
      end
 
 
 
c--------------------------------------------------------------------
c
c
      subroutine getext_i2(pix,NDIMX,NDIMY,i)
      implicit none
 
      integer NDIMX, NDIMY
      integer*2 PIX(NDIMX,NDIMY)
      integer i
 
      integer o
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer ifirst, i1, i2, ios
      integer np1, np2, npt
      integer nbu, npu
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      byte      b(2)
      integer*2 ii
      equivalence(ii,b)
 
      integer ipix, jpix
      logical islinux
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 2*NDIMX*NDIMY
      npt    =   NDIMX*NDIMY
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/2 + 1
         np2 = (nbyteE-nbyte1)/2 + 1
         read(10,rec=i,iostat=ios) buffc
         do o = 1, 1440
            npu = np1+o-1
            nbu = (o-1)*2
            ii  = 0
            ipix = npu - (npu-1)/NDIMX*NDIMX
            jpix =   1 + (npu-1)/NDIMX
            if (ipix.ge.1.and.ipix.le.NDIMX.and.
     .          jpix.ge.1.and.jpix.le.NDIMY) then
                if (.not.(islinux())) then
                   b(1) = buffb(nbu+1)
                   b(2) = buffb(nbu+2)
                   endif
                if ((islinux())) then
                   b(2) = buffb(nbu+1)
                   b(1) = buffb(nbu+2)
                   endif
                pix(ipix,jpix) = ii
                endif
            enddo!o
         enddo!i
 
      i = i2 + 1
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/xryr2mc_stdgdc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
c mc should be the pixel area correction
c
c
      subroutine xryr2mc_stdgdc(xr,yr,mc,FILE_STDGC)
      implicit none
 
      real*8 xr, yr
      real*8 mc
      character*(*) FILE_STDGC
      real*8 xc, yc
 
      call extract_stdgdc(xr,yr,xc,yc,mc,FILE_STDGC,0)  ! 0 = pix area correction
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/xryr2xcyc_stdgdc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
      subroutine xryr2xcyc_stdgdc(xr,yr,xc,yc,FILE_STDGC)
      implicit none
 
      real*8 xr, yr
      real*8 xc, yc
      character*(*) FILE_STDGC
      real*8 mc
 
      call extract_stdgdc(xr,yr,xc,yc,mc,FILE_STDGC,1)  ! 1 = forward
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/GDC/STDGDCs/xcyc2xryr_stdgdc.f"
c****                                         
c*********************************************
                                              
 
c------------------------------------------------------
c
      subroutine xcyc2xryr_stdgdc(xc,yc,xr,yr,FILE_STDGC)
      implicit none
 
      real*8 xc, yc
      real*8 xr, yr
      character*(*) FILE_STDGC
      real*8 mc
 
      call extract_stdgdc(xr,yr,xc,yc,mc,FILE_STDGC,-1)  ! -1 = reverse
 
      return
      end
 
 
 
 
 
      program dex_no_gaia
      implicit none
 
 
c
c center of field
c
      real*8  RA_CEN
      real*8  DE_CEN
      real*4  RAD
 
c
c index for GAIA stars, and data from GAIA stars
c
      integer G, Gs
      integer NREF
      character junk
      real*8      ra_g(_GMAX_)
      real*8      de_g(_GMAX_)
      real*8      er_g(_GMAX_)
      real*8      ed_g(_GMAX_)
      real*8      mg_g(_GMAX_)
      real*8     emg_g(_GMAX_)
      real*8      mb_g(_GMAX_)
      real*8     emb_g(_GMAX_)
      real*8      mr_g(_GMAX_)
      real*8     emr_g(_GMAX_)
      real*8    mura_g(_GMAX_)
      real*8   emura_g(_GMAX_)
      real*8    mude_g(_GMAX_)
      real*8   emude_g(_GMAX_)
      real*8      pi_g(_GMAX_)
      real*8     epi_g(_GMAX_)
      integer*8   id_g(_GMAX_)
      byte       b85_g(85,_GMAX_)
      real        sz_g(_GMAX_)
      integer     um_g(_GMAX_)
 
      real*8      uu_g(_GMAX_)    ! master frame coords
      real*8      vv_g(_GMAX_)    !    for gaia stars
 
      integer Os
      logical DIAG
 
      real*8 rd2x, dra            ! functions and temp
      real*8 rd2y, dde            !    variables for tan-plane work
 
      real*8 duu, dvv, ddd        ! temp variables for master frame
 
      integer P, Ps               ! info for various formulations
      integer K, Ks
      integer N, Ns
      integer L, Ls_n(4)
 
      integer Ls, Lu
      real*8 uu_l(_NSTMAX_), vv_l(_NSTMAX_), ww_l(_NSTMAX_)
      real*8  u_l(_NSTMAX_),  v_l(_NSTMAX_),  w_l(_NSTMAX_)
      real*8 ra_l(_NSTMAX_), dec_l(_NSTMAX_)
      real*8 ra_bar, dec_bar
      real*8 ra_sig, dec_sig
 
      real*8  x_ln(_NSTMAX_,4)
      real*8  y_ln(_NSTMAX_,4)
      real*8  m_ln(_NSTMAX_,4)
      real*8  p_ln(_NSTMAX_,4)
      real*8  q_ln(_NSTMAX_,4)
      real*8  u_ln(_NSTMAX_,4)
      real*8  v_ln(_NSTMAX_,4)
      real*8  w_ln(_NSTMAX_,4)
      real*8  r_ln(_NSTMAX_,4)
      real*8  d_ln(_NSTMAX_,4)
      real*8 uu_ln(_NSTMAX_,4)
      real*8 vv_ln(_NSTMAX_,4)
      real*8 ra_ln(_NSTMAX_,4)
      real*8 dec_ln(_NSTMAX_,4)
 
      real*8 quu_n(4)
      real*8 qvv_n(4)
      real*8 qmm_n(4)
 
      integer NARG, NARGs
      character*200 STRING_ARG(99)
      character*200 FILENAME
      character*200 DATADIR
      character*200 THORDIR
      character*200 GDCDIR
      character*200 PSFDIR
      integer NARG_V1, NARG_V2
      integer NARG_I1, NARG_I2
 
      logical str_contains
      character*200 STRING200
      character*999 STRING999
 
      real*8 RACEN_G, DECEN_G
      real*8 UUCEN_G, VVCEN_G
      real*8    PS_G,    TT_G
 
      integer      PID_N(4)
      character*5 FILT_N(4)
      real        EXPT_N(4)
      real*8      RDAT_N(4), RDAT
      real          ZP_N(4)
 
      real*8  BDRY_U_PK(4,4),  BDRY_V_PK(4,4)
      real*8  BDRY_R_PK(4,4),  BDRY_D_PK(4,4)
 
      real*8  RMIN, RMAX, RCEN, RWID
      real*8  DMIN, DMAX, DCEN, DWID
 
      integer hist(201,201)
      integer i, j
      integer ii,jj,hh
 
      integer U, Us, Uu
      real*8  u1_u(_NMATMX_), v1_u(_NMATMX_)
      real*8  x1_u(_NMATMX_), y1_u(_NMATMX_), m1_u(_NMATMX_)
      real*8  x2_u(_NMATMX_), y2_u(_NMATMX_), m2_u(_NMATMX_)
      real*8  x2e, y2e
      integer uu_u(_NMATMX_)
      integer  l_u(_NMATMX_)
      integer  g_u(_NMATMX_)
      real*8   ue, ve
 
      real*8 AG, BG, CG, DG
      real*8 x1o, y1o, x2o, y2o
      real   DX, DY
      real   RMS2D
 
      real*8  AG_N(4),  BG_N(4) , CG_N(4),  DG_N(4)
      real*8  GA_N(4),  GB_N(4) , GC_N(4),  GD_N(4)
      real*8 x1o_N(4), y1o_N(4), x2o_N(4), y2o_N(4)
      real   RMS_N(4)
      integer Us_N(4)
      integer Uu_N(4)
 
 
c
c HST master frame
c
      integer Q, Qs, Qu
      real*8  uu_q(_NSTMAX_),  vv_q(_NSTMAX_)
      real*8  mv_q(_NSTMAX_),  mi_q(_NSTMAX_)
      real*8 euu_q(_NSTMAX_), evv_q(_NSTMAX_)
      real*8 emv_q(_NSTMAX_), emi_q(_NSTMAX_)
      integer nv_q(_NSTMAX_)
      integer ni_q(_NSTMAX_)
      integer Lv
      integer Li
 
      real*8 uu_bar, uu_sig
      real*8 vv_bar, vv_sig
      real*8 mv_bar, mv_sig, mv_l(4)
      real*8 mi_bar, mi_sig, mi_l(4)
      real*8 ra_q(_NSTMAX_), dec_q(_NSTMAX_)
      real*8 era_q(_NSTMAX_), edec_q(_NSTMAX_)
 
 
      byte    cmap(10001,10001), c
      integer hmap(10001,10001)
      integer ih, jh
      integer ihh,jhh
      logical inside_poly
 
      integer  l_qn(_NSTMAX_,4) ! cross-link
 
      real pixi(4096,4096,2)
 
      real, dimension(:,:), allocatable :: pixo
 
      character*20 STREAM
      integer INST, INSTu
      integer V1I2
      integer NIM
      real    SCL_N(2)
      character*200 FILE_STDGDC
      integer ir, jr
      real*8  xr, yr
      real*8  xc, yc
      real*8  ug, vg
      integer ig, jg
      integer UMIN, UMAX, UDIM
      integer VMIN, VMAX, VDIM
      real    temp
      real    p_l(2)
 
      real*8  CRPIX1,CRPIX2,CRVAL1,CRVAL2
      real*8  CD1_1, CD1_2, CD2_1, CD2_2
      real*8  LTV1,  LTV2
      real*8  LTM1_1,LTM1_2,LTM2_1,LTM2_2
 
      real siz
 
c
c-----------------------------------------------------
c
 
      NARG_V1 = 0
      NARG_V2 = 0
      NARG_I1 = 0
      NARG_I2 = 0
 
      Ns = 4
 
      NARGs = iargc()
 
      if (NARGs.eq.0) then
          print*,' '
          print*,'You need to provide the routine a (potentially) long '
          print*,'   list of filenames on the command line.  It will  '
          print*,'   search through this list for things with a filter'
          print*,'   ("F606W" or "F814W") and "xympquvwrd" in the     '
          print*,'   pathname.                                        '
          print*,' '
          print*,'It will run with these files and will generate a    '
          print*,'   catalog. '
          print*,' '
          print*,'It will also search for the images associated with '
          print*,'   xympquvwrd files and will stack them, F606W and '
          print*,'   F814W.'
          print*,' '
          stop
          endif

c
c     determine the directory that holds the F606W/ and F814W/ image
c     sub-directories.  it defaults to ../DATA/ but may be overridden
c     on the command line with:  -datadir <path>
c
      DATADIR = '../DATA/'
      do NARG = 1, NARGs
         call getarg(NARG,STRING200)
         if (STRING200(1:8).eq.'-datadir') then
             if (NARG.eq.NARGs) then
                 print*,'-datadir flag given without a path.'
                 stop
                 endif
             call getarg(NARG+1,DATADIR)
             endif
         enddo!NARG
      L = len_trim(DATADIR)
      if (DATADIR(L:L).ne.'/') DATADIR = DATADIR(1:L) // '/'
      write(*,'(''        DATA DIRECTORY: '',a)') trim(DATADIR)

c
c     determine the directory that holds the geometric-distortion
c     correction (GDC) reference files.  precedence (low to high):
c       1. the relative location populated by `pip install` / setup.py,
c       2. the THOR_GDC_DIR environment variable, if set,
c       3. a -gdcdir <path> flag on the command line.
c     only the per-filter basename is appended later, so no machine-
c     specific path is ever baked into the executable.
c
      GDCDIR = '../../../src/thor/GDCs/'
      call getenv('THOR_GDC_DIR',STRING200)
      if (len_trim(STRING200).gt.0) GDCDIR = STRING200
      do NARG = 1, NARGs
         call getarg(NARG,STRING200)
         if (STRING200(1:7).eq.'-gdcdir') then
             if (NARG.eq.NARGs) then
                 print*,'-gdcdir flag given without a path.'
                 stop
                 endif
             call getarg(NARG+1,GDCDIR)
             endif
         enddo!NARG
      L = len_trim(GDCDIR)
      if (GDCDIR(L:L).ne.'/') GDCDIR = GDCDIR(1:L) // '/'
      write(*,'(''        GDC  DIRECTORY: '',a)') trim(GDCDIR)

      write(*,*)
      write(*,*)
      write(*,'(''STEP1:  DETERMINE THE PHOTOMETRY FILES'')')
      write(*,'(''        FOR EACH EXPOSURE/FILTER.     '')')
      write(*,*)
 
      do NARG = 1, NARGs
         call getarg(NARG,STRING_ARG(NARG))
         if (str_contains(STRING_ARG(NARG),200,
     .                    'F606W',5).and.
     .       str_contains(STRING_ARG(NARG),200,
     .                    'xympquvwrd',10)) then  ! it is one of the filenames!
             if (NARG_V1.eq.0) then
                 NARG_V1 = NARG
                else
                 if (NARG_V2.ne.0) then
                     print*,'You can only have two F606W files:'
                     print*,'You already have two:             '
                     write(*,'(''NARG '',i2,1x,a80)')
     .                         NARG_V1,STRING_ARG(NARG_V1)(1:80)
                     write(*,'(''NARG '',i2,1x,a80)')
     .                         NARG_V2,STRING_ARG(NARG_V2)(1:80)
                     write(*,'(''NARG '',i2,1x,a80)')
     .                         NARG, STRING_ARG(NARG)(1:80)
                     print*,' '
                     stop
                     endif
                NARG_V2 = NARG
                endif
             endif
         if (str_contains(STRING_ARG(NARG),200,
     .                    'F814W',5).and.
     .       str_contains(STRING_ARG(NARG),200,
     .                    'xympquvwrd',10)) then
             if (NARG_I1.eq.0) then
                 NARG_I1 = NARG
                else
                 if (NARG_I2.ne.0) then
                     print*,'You can only have two F814W files:'
                     print*,'You already have two:             '
                     write(*,'(''NARG '',i2,1x,a80)')
     .                         NARG_I1,STRING_ARG(NARG_I1)(1:80)
                     write(*,'(''NARG '',i2,1x,a80)')
     .                         NARG_I2,STRING_ARG(NARG_I2)(1:80)
                     write(*,'(''NARG '',i2,1x,a80)')
     .                         NARG, STRING_ARG(NARG)(1:80)
                     print*,' '
                     stop
                     endif
                NARG_I2 = NARG
                endif
             endif
         enddo!NARG
 
      write(*,*)
      write(*,'(10x,''F606W-1 NARG'',i2,1x,a80)')
     .                NARG_V1, STRING_ARG(NARG_V1)(1:80)
      write(*,'(10x,''F606W-2 NARG'',i2,1x,a80)')
     .                NARG_V2, STRING_ARG(NARG_V2)(1:80)
      write(*,*)
      write(*,'(10x,''F814W-1 NARG'',i2,1x,a80)')
     .                NARG_I1, STRING_ARG(NARG_I1)(1:80)
      write(*,'(10x,''F814W-2 NARG'',i2,1x,a80)')
     .                NARG_I2, STRING_ARG(NARG_I2)(1:80)
      write(*,*)
 
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP2:  READ IN THE STAR LISTS FROM   '')')
      write(*,'(''        THE PHOTOMETRY FILES FOR EACH '')')
      write(*,'(''        EXPOSURE.                     '')')
      write(*,*)
 
      open(21,file='dex_no_gaia_STEP02_A.nim_info',status='unknown')
      do N = 1, 4
         if (N.eq.1) FILENAME = STRING_ARG(NARG_V1)
         if (N.eq.2) FILENAME = STRING_ARG(NARG_V2)
         if (N.eq.3) FILENAME = STRING_ARG(NARG_I1)
         if (N.eq.4) FILENAME = STRING_ARG(NARG_I2)
         open(12,file=FILENAME,status='old')
         Ls_n(N) = 0
   21    continue
         read(12,'(a200)',end=22) STRING200
         if (STRING200(1:11).eq.'#      PID:')
     .       read(STRING200(12:32),*) PID_N(N)
         if (STRING200(1:11).eq.'#     FILT:')
     .       read(STRING200(12:32),*) FILT_N(N)
         if (STRING200(1:11).eq.'#    EXPTE:')
     .       read(STRING200(12:32),*) EXPT_N(N)
         if (STRING200(1:11).eq.'#     rDAT:')
     .       read(STRING200(12:32),*) RDAT_N(N)
         if (STRING200(1:13).eq.'# BDRY_uv_K1:')
     .       read(STRING200(14:200),*)  BDRY_U_PK(1,1),  BDRY_V_PK(1,1),
     .                                  BDRY_U_PK(2,1),  BDRY_V_PK(2,1),
     .                                  BDRY_U_PK(3,1),  BDRY_V_PK(3,1),
     .                                  BDRY_U_PK(4,1),  BDRY_V_PK(4,1)
         if (STRING200(1:13).eq.'# BDRY_uv_K2:')
     .       read(STRING200(14:200),*)  BDRY_U_PK(1,2),  BDRY_V_PK(1,2),
     .                                  BDRY_U_PK(2,2),  BDRY_V_PK(2,2),
     .                                  BDRY_U_PK(3,2),  BDRY_V_PK(3,2),
     .                                  BDRY_U_PK(4,2),  BDRY_V_PK(4,2)
         if (STRING200(1:13).eq.'# BDRY_rd_K1:')
     .       read(STRING200(14:200),*) BDRY_R_PK(1,1), BDRY_D_PK(1,1),
     .                                 BDRY_R_PK(2,1), BDRY_D_PK(2,1),
     .                                 BDRY_R_PK(3,1), BDRY_D_PK(3,1),
     .                                 BDRY_R_PK(4,1), BDRY_D_PK(4,1)
         if (STRING200(1:13).eq.'# BDRY_rd_K2:')
     .       read(STRING200(14:200),*) BDRY_R_PK(1,2), BDRY_D_PK(1,2),
     .                                 BDRY_R_PK(2,2), BDRY_D_PK(2,2),
     .                                 BDRY_R_PK(3,2), BDRY_D_PK(3,2),
     .                                 BDRY_R_PK(4,2), BDRY_D_PK(4,2)
         if (STRING200(1:1).eq.'#') goto 21
         Ls_n(N) = Ls_n(N)+1
         read(STRING200,*)  x_ln(Ls_N(N),N),  y_ln(Ls_N(N),N),
     .                      m_ln(Ls_N(N),N),  p_ln(Ls_N(N),N),
     .                      q_ln(Ls_N(N),N),
     .                      u_ln(Ls_N(N),N),  v_ln(Ls_N(N),N),
     .                      w_ln(Ls_N(N),N),
     .                      r_ln(Ls_N(N),N),  d_ln(Ls_N(N),N)
         goto 21
   22    continue
         close(12)
         ZP_N(N) = -2.5*log10(EXPT_N(N)/340.00)
         if (N.eq.1) write( *,221)
         if (N.eq.1) write(21,221)
         write( *,121) N, Ls_n(N),  PID_N(N), FILT_N(N),
     .                              EXPT_N(N), RDAT_N(N),
     .                              ZP_N(N),
     .                              FILENAME(1:40)
         write(21,121) N, Ls_n(N),  PID_N(N), FILT_N(N),
     .                              EXPT_N(N), RDAT_N(N),
     .                              ZP_N(N),
     .                              FILENAME(1:40)
  121    format(10x,i1,1x,i6,1x,1x,i5.5,1x,a5,1x,f6.1,1x,f12.7,
     .                 1x,f6.3,
     .                 1x,40a)
  221    format(10x,'N',1x,' Ls(N)',1x,1x,' PID ',1x,' FILT',
     .                  1x,' EXPT ',1x,'RDAT        ',
     .                  1x,' ZP   ',
     .                  1x,'FILENAME')
         enddo
      close(21)
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP3:  DETERMINE THE BOUNDARIES '')')
      write(*,*)
 
  239 format(10x,' ','N',1x,' ',i1,1x,' ',i1,
     .            5x,'  XCOR   ',
     .               '  YCOR   ',
     .            5x,'   RA       ',
     .               '   Dec      ')
      write(*,239)
      do N = 1, 4
         print*,' '
         do K = 1, 2
            do P = 1, 4
               write(*,139) N, P, K, BDRY_U_PK(P,K),
     .                               BDRY_V_PK(P,K),
     .                               BDRY_R_PK(P,K),
     .                               BDRY_D_PK(P,K)
  139          format(10x,'N',i1,1x,'P',i1,1x,'K',i1,
     .                5x,2f9.3,5x,2f12.6)
               enddo
            enddo
         enddo
      print*,'     '
 
      RMIN =  999.99
      RMAX = -999.99
      DMIN =  999.99
      DMAX = -999.99
      do N = 1, Ns
         do L = 1, Ls_N(N)
            if (r_ln(L,N).lt.RMIN) RMIN = r_ln(L,N)
            if (r_ln(L,N).gt.RMAX) RMAX = r_ln(L,N)
            if (d_ln(L,N).lt.DMIN) DMIN = d_ln(L,N)
            if (d_ln(L,N).gt.DMAX) DMAX = d_ln(L,N)
            enddo
         enddo
 
      RCEN = (RMIN+RMAX)/2
      RWID = (RMAX-RMIN)
 
      DCEN = (DMIN+DMAX)/2
      DWID = (DMAX-DMIN)
 
      write(*,*)                      !1234567890123
  449 format(10x,'     ',a13,a13,a13,a13)
      write(*,449) ' ............',' ............',
     .             ' ............',' ............'
      write(*,449) '     MIN     ','     MAX     ',
     .             '     CEN     ','     WID     '
      write(*,449) ' ............',' ............',
     .             ' ............',' ............'
      write(*,'(10x,'' RA: '',4f13.6)') RMIN,RMAX,RCEN,RWID
      write(*,'(10x,''Dec: '',4f13.6)') DMIN,DMAX,DCEN,DWID
      write(*,*)
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP4:  INTERNAL REF STARS (NO GAIA DATABASE)'')')
      write(*,'(''       Bright stars from exposure NREF'')')
      write(*,*)
 
      RA_CEN = RCEN
      DE_CEN = DCEN
 
      NREF = 1
 
      Gs = 0
      do L = 1, Ls_n(NREF)
         if (m_ln(L,NREF).lt.-11) then
             Gs = Gs + 1
             if (Gs.gt._GMAX_) stop 'dex_no_gaia: Gs.gt._GMAX_'
             ra_g(Gs) = r_ln(L,NREF)
             de_g(Gs) = d_ln(L,NREF)
             id_g(Gs) = int(L, kind=8)
             mg_g(Gs) = 18.d0
             mb_g(Gs) = 18.d0
             mr_g(Gs) = 18.d0
             emg_g(Gs) = 0.01d0
             emb_g(Gs) = 0.01d0
             emr_g(Gs) = 0.01d0
             er_g(Gs) = 0.01d0
             ed_g(Gs) = 0.01d0
             mura_g(Gs) = 0.d0
             mude_g(Gs) = 0.d0
             emura_g(Gs) = 1.d0
             emude_g(Gs) = 1.d0
             pi_g(Gs) = 0.d0
             epi_g(Gs) = 1.d0
             endif
      enddo
 
      write(*,*)
      write(*,'(10x,''    Gs = '',i6)') Gs
      write(*,*)
      print*,'     '
 
      RDAT = RDAT_N(1)
      print*,'  RDAT = ',RDAT
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP5:  TRANSFORM INTO REF FRAME  '')')
      write(*,*)
      write(*,'(10x,''  RDAT = '',f8.2)') RDAT
      write(*,*)
 
 
      open(51,file='dex_no_gaia_STEP05_A.uvm_gaia',status='unknown')
      write(51,'(''#'')')
      write(51,'(''# RA_CEN: '',f15.8)') RA_CEN
      write(51,'(''# DE_CEN: '',f15.8)') DE_CEN
      write(51,'(''#'')')
      write(51,'(''#  COL(01): ugaia - gaia map to ref frame (x) '')')
      write(51,'(''#  COL(02): vgaia - gaia map to ref frame (y) '')')
      write(51,'(''#  COL(03): mgaia - gaia gmag '')')
      write(51,'(''#  COL(04): bgaia - gaia bmag '')')
      write(51,'(''#  COL(05): rgaia - gaia rmag '')')
      write(51,'(''#  COL(06): ra0_g - gaia RA (nopm)       '')')
      write(51,'(''#  COL(07): de0_g - gaia Dec (nopm)      '')')
      write(51,'(''#  COL(08): unopm - gaia map w/o pm (x)  '')')
      write(51,'(''#  COL(09): vnopm - gaia map w/o pm (y)  '')')
      write(51,'(''#'')')
      do G = 1, Gs
         dra = rd2x(ra_g(G),de_G(G),RA_CEN,DE_CEN)
         dde = rd2y(ra_g(G),de_G(G),RA_CEN,DE_CEN)
c        uu_g(G) = 5001 - dra*60*60/0.040                               ! without any correction for proper
c        vv_g(G) = 5001 + dde*60*60/0.040                               !     motion
         uu_g(G) = 5001 - dra*60*60/0.040 - (RDAT-2016.0)*mura_g(G)/40  ! with correction for proper
         vv_g(G) = 5001 + dde*60*60/0.040 + (RDAT-2016.0)*mude_g(G)/40  !    motion (mura/mude are in mas/yr)
         if (mg_g(G).gt.1)
     .   write(51,151) uu_g(G), vv_g(G), mg_g(G), mg_g(G), mr_g(G),
     .                 ra_g(G), de_g(G),
     .                 5001 - dra*60*60/0.040,                          ! without any correction for proper
     .                 5001 + dde*60*60/0.040                           !     motion
  151    format(1x,f9.3,1x,f9.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,
     .          1x,f15.8,1x,f15.8,1x,
     .          1x, f9.3,1x, f9.3)
         enddo
      close(51)
 
      open(52,file='dex_no_gaia_STEP05_B.uvrd_1p',status='unknown')
      write(52,'(''#'')')
      write(52,'(''# RA_CEN: '',f15.8)') RA_CEN
      write(52,'(''# DE_CEN: '',f15.8)') DE_CEN
      write(52,'(''#'')')
      do N = 1, 4
         do L = 1, Ls_N(N)
            dra = rd2x(r_ln(L,N),d_ln(L,N),RA_CEN,DE_CEN)
            dde = rd2y(r_ln(L,N),d_ln(L,N),RA_CEN,DE_CEN)
            uu_ln(L,N) = 5001 - dra*60*60/0.040
            vv_ln(L,N) = 5001 + dde*60*60/0.040
            write(52,152) uu_ln(L,N), vv_ln(L,N),
     .                     r_ln(L,N),  d_ln(L,N), L, N
  152       format(1x,f9.4,1x,f9.4,1x,f15.8,1x,f15.8,' L ',i5,' N ',i1)
            enddo
         write(*,'(10x,''   N'',i1,''s = '',i6)') N,Ls_N(N)
         enddo
      close(52)
 
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP6:  MATCH EACH EXPOSURE TO INTERNAL REF'')')
      write(*,*)
 
      open(61,file='dex_no_gaia_STEP06_A.mat_gaia_list',
     .        status='unknown')
      open(62,file='dex_no_gaia_STEP06_B.mat_gaia_info',
     .        status='unknown')
      write(61,'(''#'')')
      write(61,'(''#  COL(01):  xgc  '')')
      write(61,'(''#  COL(02):  ygc  '')')
      write(61,'(''#  COL(03): xgaia '')')
      write(61,'(''#  COL(04): ygaia '')')
      write(61,'(''#  COL(05):  mgc  '')')
      write(61,'(''#  COL(06): mgaia '')')
      write(61,'(''#  COL(07):   DX  '')')
      write(61,'(''#  COL(08):   DY  '')')
      write(61,'(''#  COL(09):  x2e  '')')
      write(61,'(''#  COL(10):  y2e  '')')
      write(61,'(''#  COL(11):   N   '')')
      write(61,'(''#'')')
      write(62,262)
      write(62,'(''#'')')
      print*,' '
      do N = 1, 4
         Ls = 0
         do L = 1, Ls_N(N)
            if (m_ln(L,N).lt.-11) then
                Ls = Ls + 1
                uu_l(Ls) = uu_ln(L,N)
                vv_l(Ls) = vv_ln(L,N)
                ww_l(Ls) =  w_ln(L,N)
                 u_l(Ls) =  u_ln(L,N)
                 v_l(Ls) =  v_ln(L,N)
                 w_l(Ls) =  w_ln(L,N)
                endif
            enddo
         write(*,*)
         write(*,'(10x,''   N: '',i1,1x,i6,1x,i6)')
     .                             N, Ls_N(N), Ls
         write(*,*)
         do i = -100, 100
         do j = -100, 100
            hist(101+i,101+j) = 0
            enddo
            enddo
         do L = 1, Ls
         do G = 1, Gs
            duu = uu_l(L)-uu_g(G)
            dvv = vv_l(L)-vv_g(G)
            i = 101 + duu + 0.5
            j = 101 + dvv + 0.5
            if (i.ge.002.and.i.le.200.and.
     .          j.ge.002.and.j.le.200) then
                hist(i  ,j  ) = hist(i  ,j  ) + 1
                hist(i+1,j  ) = hist(i+1,j  ) + 1
                hist(i-1,j  ) = hist(i-1,j  ) + 1
                hist(i  ,j+1) = hist(i  ,j+1) + 1
                hist(i  ,j-1) = hist(i  ,j-1) + 1
                hist(i+1,j+1) = hist(i+1,j+1) + 1
                hist(i+1,j-1) = hist(i+1,j-1) + 1
                hist(i-1,j+1) = hist(i-1,j+1) + 1
                hist(i-1,j-1) = hist(i-1,j-1) + 1
                endif
            enddo
            enddo
         write(*,*)
         do j = 10, -10, -1
            write(*,'(10x,i3,1x,21i5)') j, (hist(101+i,101+j),i=-10,10)
            enddo
         write(*,*)
         ii = 0
         jj = 0
         hh = 0
         do i = 001, 201
         do j = 001, 201
            if (hist(i,j).gt.hh) then
                hh = hist(i,j)
                ii = i
                jj = j
                endif
            enddo
            enddo
         write(*,*)
         write(*,'(10x,''   HMAX: '',i5,1x,i3,1x,i3)')
         write(*,*)
         Us = 0
         do L = 1, Ls
         do G = 1, Gs
            duu = uu_l(L)-uu_g(G)
            dvv = vv_l(L)-vv_g(G)
            i = 101 + duu + 0.5
            j = 101 + dvv + 0.5
            if ((i-ii)**2+(j-jj)**2.lt.2.5**2.and.
     .           mg_g(G).gt.01.and.
     .           mg_g(G).lt.19) then
                 Us = Us + 1
                 if (Us.gt.99999) stop 'Us.gt.9999'
                 u1_u(Us) = uu_l(L)
                 v1_u(Us) = vv_l(L)
                 x1_u(Us) =  u_l(L)
                 y1_u(Us) =  v_l(L)
                 m1_u(Us) =  w_l(L)
                 x2_u(Us) = uu_g(G)
                 y2_u(Us) = vv_g(G)
                 m2_u(Us) = mg_g(G)
                 uu_u(Us) = 1
                  l_u(Us) = L
                  g_u(Us) = G
                 endif
            enddo
            enddo
         write(*,*)
         write(*,'(10x,''   HMAX: '',i5,1x,i5,1x,i5)')
     .                               Us,   Ls,   Gs
         call glob_fit6nrDPU_REJ(u1_u,v1_u,x2_u,y2_u,uu_u,Us,Uu,           ! added Nu
     .                           AG,BG,CG,DG,x1o,y1o,x2o,y2o,
     .                           RMS2D,2.75,0.05,1.00) ! added sev
         DX = x1o-x2o
         DY = y1o-y2o
         write(*,*)
         write(*,'(10x,''     Uu: '',i7,'' / '', i7)') Uu, Us
         write(*,'(10x,''     AB: '',f12.7,1x,f12.7)') AG, BG
         write(*,'(10x,''     CD: '',f12.7,1x,f12.7)') CG, DG
         write(*,'(10x,''     DX: '',f12.3,1x,f12.3)') DX, DY
         write(*,'(10x,''    RMS: '',f12.7)') RMS2D
 
         write(*,*)
         call glob_fit6nrDPU_REJ(x1_u,y1_u,x2_u,y2_u,uu_u,Us,Uu,           ! added Nu
     .                           AG,BG,CG,DG,x1o,y1o,x2o,y2o,
     .                           RMS2D,2.75,0.05,1.00) ! added sev
         write(*,*)
         write(*,'(10x,''     Uu: '',i7,'' / '', i7)') Uu, Us
         write(*,'(10x,''     AB: '',f12.7,1x,f12.7)') AG, BG
         write(*,'(10x,''     CD: '',f12.7,1x,f12.7)') CG, DG
         write(*,'(10x,''     DX: '',f12.3,1x,f12.3)') DX, DY
         write(*,'(10x,''    RMS: '',f12.7)') RMS2D
 
         AG_N(N) = AG
         BG_N(N) = BG
         CG_N(N) = CG
         DG_N(N) = DG
         GA_N(N) =  DG/(AG*DG-BG*CG)
         GB_N(N) = -BG/(AG*DG-BG*CG)
         GC_N(N) = -CG/(AG*DG-BG*CG)
         GD_N(N) =  AG/(AG*DG-BG*CG)
         x1o_N(N) = x1o
         y1o_N(N) = y1o
         x2o_N(N) = x2o
         y2o_N(N) = y2o
         RMS_N(N) = RMS2D
         Us_N(N) = Us
         Uu_N(N) = Uu
         do U = 1, Us
            if (uu_u(U).eq.1) then
               x2e = x2o + AG*(x1_u(U)-x1o) + BG*(y1_u(U)-y1o)
               y2e = y2o + CG*(x1_u(U)-x1o) + DG*(y1_u(U)-y1o)
               write(61,161) x1_u(U), y1_u(U), x2_u(U), y2_u(U),
     .                       m1_u(U), m2_u(U),
     .                       x2_u(U)-x2e, y2_u(U)-y2e, x2e, y2e, N
  161          format(1x,f9.3,1x,f9.3,3x,f9.3,1x,f9.3,
     .                3x,f8.3,1x,f8.3,
     .                3x,f8.3,1x,f8.3,
     .                5x,f9.3,1x,f9.3,
     .                1x,i1,' N')
               endif
            enddo
         enddo
 
      print*,' '
  124 format('#','N',5x,'.Us.',1x,'.Uu.',5x,
     .           ' ..AG......',
     .           ' ..BG......',
     .           ' ..CG......',
     .           ' ..DG......',5x,
     .           ' ..RMS..')
      do N = 1, 4
         write(62,162) N, Us_N(N), Uu_N(N),
     .                 AG_N(N), BG_N(N), CG_N(N), DG_N(N),
     .                 x1o_N(N)+GA_N(N)*(5001.0-x2o_N(N))
     .                         +GB_N(N)*(5001.0-y2o_N(N)),
     .                 y1o_N(N)+GC_N(N)*(5001.0-x2o_N(N))
     .                         +GD_N(N)*(5001.0-y2o_N(N)),
     .                 5001.0, 5001.0,
     .                 RMS_N(N), DX, DY
         write( *,162) N, Us_N(N), Uu_N(N),
     .                 AG_N(N), BG_N(N), CG_N(N), DG_N(N),
     .                 x1o_N(N)+GA_N(N)*(5001.0-x2o_N(N))
     .                         +GB_N(N)*(5001.0-y2o_N(N)),
     .                 y1o_N(N)+GC_N(N)*(5001.0-x2o_N(N))
     .                         +GD_N(N)*(5001.0-y2o_N(N)),
     .                 5001.0, 5001.0,
     .                 RMS_N(N), DX, DY
  162    format(2x,i1,5x,i5,1x,i5,
     .                5x,4f11.7,5x,2f9.3,
     .                5x,2f9.3,
     .                5x,f8.4,1x,f7.3,1x,f7.3)
  262    format('# ','N',5x,'  Us ',1x,'  Uu ',
     .                   5x,'    AG     ','    BG     ',
     .                      '    CG     ','    DG     ',
     .                   5x,'  xgc0   ','  ygc0   '
     .                   5x,'   u0    ','   v0    '
     .                   5x,'  RMS   ',1x,'   DX  ',1x,'   DX  ')
         enddo
      print*,' '
      close(61)
      close(62)
 
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP7:  MAKE A MASTER LIST OF ALL THE STARS'')')
      write(*,*)
 
c----------------------------------------------------------
c
c now, go through and make a master list of all the stars
c
      do ih = 0001, 10001
      do jh = 0001, 10001
         hmap(ih,jh) = 0
         cmap(ih,jh) = 0
         enddo
         enddo
 
      do N = 1, Ns
         do ih = 0001, 10001
         do jh = 0001, 10001
            ue = x1o_N(N) + GA_N(N)*(ih-x2o_n(N))
     .                    + GB_N(N)*(jh-y2o_n(N))
            ve = y1o_N(N) + GC_N(N)*(ih-x2o_n(N))
     .                    + GD_N(N)*(jh-y2o_n(N))
            if (inside_poly(ue,ve,BDRY_U_PK(1,1),BDRY_V_PK(1,1),4))
     .          cmap(ih,jh) = cmap(ih,jh) + 1
            if (inside_poly(ue,ve,BDRY_U_PK(1,2),BDRY_V_PK(1,2),4))
     .          cmap(ih,jh) = cmap(ih,jh) + 1
            enddo
            enddo
         enddo
 
      Qs = 0
      do N = 1, Ns
         do L = 1, Ls_n(N)
            uu_ln(L,N) = x2o_N(N) + AG_N(N)*(u_ln(L,N)-x1o_N(N))
     .                            + BG_N(N)*(v_ln(L,N)-y1o_N(N))
            vv_ln(L,N) = y2o_N(N) + CG_N(N)*(u_ln(L,N)-x1o_N(N))
     .                            + DG_N(N)*(v_ln(L,N)-y1o_N(N))
            ih = int(uu_ln(L,N))
            jh = int(vv_ln(L,N))
            if (ih.ge.00001.and.ih.le.10001.and.
     .          jh.ge.00001.and.jh.le.10001) then
                hmap(ih  ,jh  ) = hmap(ih  ,jh  ) + 1
                hmap(ih+1,jh  ) = hmap(ih+1,jh  ) + 1
                hmap(ih  ,jh+1) = hmap(ih  ,jh+1) + 1
                hmap(ih+1,jh+1) = hmap(ih+1,jh+1) + 1
                endif
            enddo
         enddo!N
 
      Qs = 0
      do ih = 00001+2, 10001-2
      do jh = 00001+2, 10001-2
         if (hmap(ih,jh).ge.2.and.                  ! found in at least 3/4 exposures
     .       hmap(ih,jh).ge.hmap(ih-1,jh  ).and.
     .       hmap(ih,jh).ge.hmap(ih+1,jh  ).and.
     .       hmap(ih,jh).ge.hmap(ih  ,jh-1).and.
     .       hmap(ih,jh).ge.hmap(ih  ,jh+1).and.
     .       hmap(ih,jh).ge.hmap(ih-1,jh+1).and.
     .       hmap(ih,jh).ge.hmap(ih+1,jh+1).and.
     .       hmap(ih,jh).ge.hmap(ih-1,jh-1).and.
     .       hmap(ih,jh).ge.hmap(ih+1,jh-1)) then
             Qs = Qs + 1
             uu_q(Qs) = ih
             vv_q(Qs) = jh
             do ihh = ih-1, ih+1
             do jhh = jh-1, jh+1
                hmap(ihh,jhh) = -Qs
                enddo
                enddo
             do N = 1, 4
                l_qn(Qs,N) = 0
                enddo
             endif
         enddo
         enddo
 
      write(*,*)
      write(*,'(10x,''  MASTER LIST IN 3+/4 --> Qs: '',i6)') Qs
      write(*,*)
 
 
      write(*,*)
      write(*,*)
      write(*,'(''STEP8:  REFINE MASTER LIST '')')
      write(*,*)
 
      open(81,file='dex_no_gaia_STEP08_A.xyvieeee',
     .        status='unknown')
      do N = 1, 4
         Lu = 0
         do L = 1, Ls_n(N)
            ih = int(uu_ln(L,N)+0.5)
            jh = int(vv_ln(L,N)+0.5)
            if (hmap(ih,jh).lt.0) then
                Q = -hmap(ih,jh)
                l_qn(Q,N) = L
                Lu = Lu + 1
                endif
            enddo!L
         print*,'                         N: ',N,Lu
         enddo
      print*,' '
 
      do Q = 1, Qs
         Ls = 0
         Lv = 0
         Li = 0
         do N = 1, Ns
            quu_n(N) = 0.0
            qvv_n(N) = 0.0
            qmm_n(N) = 0.0
            if (l_qn(Q,N).ne.0) then
                Ls = Ls + 1
                uu_l(Ls) = uu_ln(l_qn(Q,N),N)
                vv_l(Ls) = vv_ln(l_qn(Q,N),N)
                ra_l(Ls)  = r_ln(l_qn(Q,N),N)
                dec_l(Ls) = d_ln(l_qn(Q,N),N)
                quu_n(N) = uu_l(Ls)
                qvv_n(N) = vv_l(Ls)
                qmm_n(N) = w_ln(l_qn(Q,N),N) - ZP_N(N)
                if (FILT_N(N).eq.'F606W') then
                    Lv = Lv + 1
                    mv_l(Lv) = w_ln(l_qn(Q,N),N) - ZP_N(N)
                    endif
                if (FILT_N(N).eq.'F814W') then
                    Li = Li + 1
                    mi_l(Li) = w_ln(l_qn(Q,N),N) - ZP_N(N)
                    endif
                endif
            enddo
c        write(99,199) Q, (quu_n(N),N=1,4),
c    .                    (qvv_n(N),N=1,4),
c    .                    (qmm_n(N),N=1,4)
c 199    format(1x,i5.5,1x,4f9.3,1x,4f9.3,1x,4f8.3)
         call dbarsigs_2d(uu_l,vv_l,Ls,
     .                    uu_bar,uu_sig,
     .                    vv_bar,vv_sig,Lu,3.25)
         mv_bar = 0.00
         mv_sig = 9.99
         if (Lv.eq.2) then
             mv_bar = (mv_l(1)+mv_l(2))/2.0
             mv_sig = abs(mv_l(1)-mv_l(2))/2.0
             endif
         if (Lv.eq.1) mv_bar = mv_l(1)
         mi_bar = 0.00
         mi_sig = 9.99
         if (Li.eq.2) then
             mi_bar = (mi_l(1)+mi_l(2))/2.0
             mi_sig = abs(mi_l(1)-mi_l(2))/2.0
             endif
         if (Li.eq.1) mi_bar = mi_l(1)
         ra_bar = 0.d0
         dec_bar = 0.d0
         ra_sig = 9.99d0
         dec_sig = 9.99d0

         if (Ls.ge.1) then
            do i = 1, Ls
               ra_bar  = ra_bar  + ra_l(i)
               dec_bar = dec_bar + dec_l(i)
            enddo
               ra_bar  = ra_bar  / dble(Ls)
               dec_bar = dec_bar / dble(Ls)
         endif

         if (Ls.ge.2) then
            ra_sig = 0.d0
            dec_sig = 0.d0
            do i = 1, Ls
                  ra_sig  = ra_sig  + (ra_l(i)  - ra_bar )**2
                  dec_sig = dec_sig + (dec_l(i) - dec_bar)**2
            enddo
            ra_sig  = sqrt(ra_sig  / dble(Ls-1))
            dec_sig = sqrt(dec_sig / dble(Ls-1))
         endif
         if (Q.eq.Q/1000*1000)
     .   write( *,181) uu_bar, vv_bar, mv_bar, mi_bar,
     .                 uu_sig, vv_sig, mv_sig, mi_sig,
     .                 ra_bar, dec_bar,
     .                 Ls, Lu, Lv, Li, Q
         write(81,181) uu_bar, vv_bar, mv_bar, mi_bar,
     .                 uu_sig, vv_sig, mv_sig, mi_sig,
     .                 ra_bar, dec_bar,
     .                 Ls, Lu, Lv, Li, Q
  181 format(
     .        1x,f8.3, 1x,f8.3, 1x,f8.3, 1x,f8.3,
     .        5x,f6.3, 1x,f6.3, 1x,f6.3, 1x,f6.3,
     .        2x,f12.7, 1x,f12.7,
     .        2x,i1, 1x,i1, 1x,i1, 1x,i1, 1x,i6.6
     .       )
         uu_q(Q) = uu_bar
         vv_q(Q) = vv_bar
         mv_q(Q) = mv_bar
         mi_q(Q) = mi_bar
         nv_q(Q) = Lv
         ni_q(Q) = Li
         euu_q(Q) = uu_sig
         evv_q(Q) = vv_sig
         emv_q(Q) = mv_sig
         emi_q(Q) = mi_sig
         ra_q(Q) = ra_bar
         dec_q(Q) = dec_bar
         era_q(Q) = ra_sig
         edec_q(Q) = dec_sig
         enddo
 
 
      print*,' '
      print*,' '
      print*,'MAT FILES... '
      print*,' '
      print*,' '
      do N = 1, Ns
         Us = 0
         do Q = 1, Qs
            L = l_qn(Q,N)
            if (L.ge.0) then
                if (q_ln(L,N).gt.  0.001.and.
     .              q_ln(L,N).lt.  0.500.and.
     .              m_ln(L,N).lt. -9.000.and.
     .              m_ln(L,N).gt.-13.750.and.
     .              nv_q(Q).eq.2.and.
     .              ni_q(Q).eq.2) then
                    Us = Us + 1
                    x1_u(Us) =  u_ln(L,N)
                    y1_u(Us) =  v_ln(L,N)
                    m1_u(Us) =  w_ln(L,N)
                    x2_u(Us) = uu_q(Q)
                    y2_u(Us) = vv_q(Q)
                    if (FILT_N(N).eq.'F606W') m2_u(Us) = mv_q(Q)
                    if (FILT_N(N).eq.'F814W') m2_u(Us) = mi_q(Q)
                    uu_u(Us) = 1
                    endif
                endif
            enddo
         call glob_fit6nrDPU_REJ(x1_u,y1_u,x2_u,y2_u,uu_u,Us,Uu,           ! added Nu
     .                           AG,BG,CG,DG,x1o,y1o,x2o,y2o,
     .                           RMS2D,2.75,0.05,1.00) ! added sev
         AG_N(N) = AG
         BG_N(N) = BG
         CG_N(N) = CG
         DG_N(N) = DG
         GA_N(N) =  DG/(AG*DG-BG*CG)
         GB_N(N) = -BG/(AG*DG-BG*CG)
         GC_N(N) = -CG/(AG*DG-BG*CG)
         GA_N(N) =  AG/(AG*DG-BG*CG)
         x1o_N(N) = x1o
         y1o_N(N) = y1o
         x2o_N(N) = x2o
         y2o_N(N) = y2o
         RMS_N(N) = RMS2D
         Us_N(N) = Us
         Uu_N(N) = Uu
         write(*,162) N, Us_N(N), Uu_N(N),
     .                AG_N(N), BG_N(N), CG_N(N), DG_N(N),
     .                x1o_N(N)+GA_N(N)*(5001.0-x2o_N(N))
     .                        +GB_N(N)*(5001.0-y2o_N(N)),
     .                y1o_N(N)+GC_N(N)*(5001.0-x2o_N(N))
     .                        +GD_N(N)*(5001.0-y2o_N(N)),
     .                 5001.0, 5001.0,
     .                RMS_N(N)
         enddo!N
 
 
      do i = 00001, 10001
      do j = 00001, 10001
         hmap(i,j) = 0
         enddo
         enddo
 
      open(91,file='dex_no_gaia_STEP09_A.gaia_linkup',
     .        status='unknown')
      write(91,391)
      write(91,291)
      write(91,391)
  391 format('#','.....',1x,'.........',1x,'.........',1x,
     .                   1x,'.......',1x,'.......',1x,'.......',
     .                   1x,'....................',1x,'.',
     .        5x,'.....',1x,'.........',1x,'.........',
     .                   1x,'.......',1x,'.......',
     .                   1x,'......',1x,'......',
     .                   1x,'......',1x,'......',1x,'.',1x,'.')
  291 format('#','    G',1x,'  U_GAIA ',1x,'  V_GAIA ',1x,
     .                   1x,'MG_GAIA',1x,'MB_GAIA',1x,'MR_GAIA',
     .                   1x,'    ID_GAIA         ',1x,'C',
     .        5x,'    Q',1x,'    U_HST',1x,'    V_HST',
     .                   1x,' MV_HST',1x,' MI_HST',
     .                   1x,'eU_HST',1x,'eV_HST',
     .                   1x,'   eMV',1x,'   eMI',1x,'V',1x,'I')
 
 
      do G = 1, Gs
         Q = 0                                     ! no star found
         C  = 0
         uu_bar = 0.
         vv_bar = 0.
         mv_bar = 0.
         mi_bar = 0.
         ra_bar = 0.d0
         dec_bar = 0.d0
         uu_sig = 9.99
         vv_sig = 9.99
         mv_sig = 9.99
         mi_sig = 9.99
         ra_sig = 9.99d0
         dec_sig = 9.99d0
         Lv = 0
         i = int(uu_g(G)+0.5)
         j = int(vv_g(G)+0.5)
         if (i.ge.00001+5.and.i.le.10001-5.and.
     .       j.ge.00001+5.and.j.le.10001-5.and.
     .       mg_g(G).gt.1.00) then
             C = cmap(i,j)
             if (C.ge.0) then              ! this GAIA star could in principle be found in at least one image
                 do Qu = 1, Qs
                    duu = uu_g(G)-uu_q(Qu)
                    dvv = vv_g(G)-vv_q(Qu)
                    ddd = sqrt(duu**2+dvv**2)
                    if (ddd.lt.0.75) then
                        Q = Qu
                        uu_bar =  uu_q(Q)
                        vv_bar =  vv_q(Q)
                        mv_bar =  mv_q(Q)
                        mi_bar =  mi_q(Q)
                        uu_sig = euu_q(Q)
                        vv_sig = evv_q(Q)
                        mv_sig = emv_q(Q)
                        mi_sig = emi_q(Q)
                        Lv     = nv_q(Q)
                        Li     = ni_q(Q)
                        ra_bar  = ra_q(Q)
                        dec_bar = dec_q(Q)
                        ra_sig  = era_q(Q)
                        dec_sig = edec_q(Q)
                        endif
                    enddo
                 endif
             write(91,191) G, uu_g(G), vv_g(G),
     .                        mg_g(G), mb_g(G), mr_g(G),
     .                        id_g(G), C,
     .                     Q, uu_bar, vv_bar, mv_bar, mi_bar,
     .                        uu_sig, vv_sig, mv_sig, mi_sig,
     .                        Lv, Li
  191        format(1x,i5,1x,f9.3,1x,f9.3,1x,
     .                    1x,f7.3,1x,f7.3,1x,f7.3,1x,i20,1x,i1,
     .              5x,i5.5,1x,f9.3,1x,f9.3,1x,f7.3,1x,f7.3,
     .                      1x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3,1x,i1,1x,i1)
             endif
         enddo
      write(91,391)
      write(91,291)
      write(91,391)
      close(91)
 
 
      UMIN = 99999
      UMAX = -9999
      VMIN = 99999
      VMAX = -9999
 
      INSTu = 0
      do V1I2 = 1, 2
         print*,' '
         if (V1I2.eq.1) print*,'READ IN IMAGES FOR F606W... '
         if (V1I2.eq.2) print*,'READ IN IMAGES FOR F814W... '
         print*,' '
         do NIM = 1, 2
            if (V1I2.eq.1.and.NIM.eq.1) FILENAME = trim(DATADIR)
     .                               // 'F606W/'
     .                               // trim(STRING_ARG(NARG_V1))
            if (V1I2.eq.1.and.NIM.eq.2) FILENAME = trim(DATADIR)
     .                               // 'F606W/'
     .                               // trim(STRING_ARG(NARG_V2))
            if (V1I2.eq.2.and.NIM.eq.1) FILENAME = trim(DATADIR)
     .                               // 'F814W/'
     .                               // trim(STRING_ARG(NARG_I1))
            if (V1I2.eq.2.and.NIM.eq.2) FILENAME = trim(DATADIR)
     .                               // 'F814W/'
     .                               // trim(STRING_ARG(NARG_I2))
            do i = 1, 180
               if (FILENAME(i:i+4).eq.'_flc.') then
                   FILENAME = FILENAME(1:i+4) // 'fits'
                   goto 108
                   endif
               enddo
            print*,'No "_flc." found in FILENAME.'
            stop   'cannot read in the images.'
  108       continue
            call readfits_r4e(FILENAME,pixi(0001,0001,NIM),4096,2048,1)
            call readfits_r4e(FILENAME,pixi(0001,2049,NIM),4096,2048,4)
            call query_hdr(FILENAME,'DETECTOR',STREAM)
            INST = 0
            if (STREAM(2:4).eq.'WFC' ) INST = 1
            if (STREAM(2:5).eq.'UVIS') INST = 2
            if (INSTu.ne.0.and.INSTu.ne.INST)
     .          stop 'inconsistent instruments!'
            INSTu = INST
            N = NIM + (V1I2-1)*2
            SCL_N(NIM) = EXPT_N(N)/340.00
            FILE_STDGDC = 'NONE'
            if (V1I2.eq.1.and.INST.eq.1)
     .          FILE_STDGDC = trim(GDCDIR)
     .                     // 'STDGDC_OFFICIAL_JFRAME_ACSWFC_F606W.fits'
            if (V1I2.eq.2.and.INST.eq.1)
     .          FILE_STDGDC = trim(GDCDIR)
     .                     // 'STDGDC_OFFICIAL_JFRAME_ACSWFC_F814W.fits'
            if (V1I2.eq.1.and.INST.eq.2)
     .          FILE_STDGDC = trim(GDCDIR)
     .                     // 'STDGDC_WFC3UV_F606W.fits'
            if (V1I2.eq.2.and.INST.eq.2)
     .          FILE_STDGDC = trim(GDCDIR)
     .                     // 'STDGDC_WFC3UV_F814W.fits'
            write(*,117) V1I2,NIM,N,INST,INSTu,
     .                   int(EXPT_N(N)+0.5),SCL_N(NIM),
     .                   pixi(2048,1024,NIM),
     .                   FILENAME(1:30),
     .                   FILE_STDGDC(1:82)
  117       format(5x,i1,1x,i1,1x,i1,1x,i1,1x,i1,
     .                   1x,i3,1x,f8.4,
     .                   1x,f8.2,
     .                   1x,a30,
     .                   1x,a82)
            if (V1I2.eq.1) then
                do i = 0001, 4096, 4095
                do j = 0001, 4096, 4095
                   xr = i
                   yr = j
                   call xryr2xcyc_STDGDC(xr,yr,xc,yc,FILE_STDGDC)
                   ug = x2o_N(N) + AG_N(N)*(xc-x1o_N(N))
     .                           + BG_N(N)*(yc-y1o_N(N))
                   vg = y2o_N(N) + CG_N(N)*(xc-x1o_N(N))
     .                           + DG_N(N)*(yc-y1o_N(N))
                   if (ug.lt.UMIN) UMIN = ug
                   if (ug.gt.UMAX) UMAX = ug
                   if (vg.lt.VMIN) VMIN = vg
                   if (vg.gt.VMAX) VMAX = vg
                   enddo
                   enddo
               LTV1   = -UMIN + 1
               LTV2   = -VMIN + 1
               LTM1_1 =  1.00
               LTM1_2 =  0.00
               LTM2_1 =  0.00
               LTM2_2 =  1.00
               CRPIX1 = 5001.00-UMIN
               CRPIX2 = 5001.00-VMIN
               CRVAL1 = RA_CEN
               CRVAL2 = DE_CEN
               CD1_1  = -1.388888904E-05*(40./50)
               CD1_2  =  0.000000000E-05*(40./50)
               CD2_1  =  0.000000000E-05*(40./50)
               CD2_2  =  1.388888904E-05*(40./50)
               endif
            enddo
         if (V1I2.eq.1) then
             UDIM = UMAX-UMIN+1
             VDIM = VMAX-VMIN+1
             allocate(pixo(UDIM,VDIM))
             endif
         print*,' '
         print*,'  U: ',UMIN,UMAX,UDIM
         print*,'  V: ',VMIN,VMAX,VDIM
         print*,' '
         do j = 1, VDIM
         do i = 1, UDIM
            pixo(i,j) = -100
            Ls = 0
            Lu = 0
            do NIM = 1, 2
               N = NIM + (V1I2-1)*2
               ig = UMIN + (i-1)
               jg = VMIN + (j-1)
               ug = ig
               vg = jg
               xc = x1o_N(N) + GA_N(N)*(ug-x2o_N(N))
     .                       + GB_N(N)*(vg-y2o_N(N))
               yc = y1o_N(N) + GC_N(N)*(ug-x2o_N(N))
     .                       + GD_N(N)*(vg-y2o_N(N))
               call xcyc2xryr_STDGDC(xc,yc,xr,yr,FILE_STDGDC)
               ir = xr+0.5
               jr = yr+0.5
               if (ir.ge.1.and.ir.le.4096.and.
     .             jr.ge.1.and.jr.le.4096.and.
     .             abs(jr-2048.5).gt.4.5) then
                   Ls = Ls + 1
                   p_l(Ls) = pixi(ir,jr,NIM)/SCL_N(NIM)
                   endif
               enddo!NIM
            if (Ls.eq.0) goto 9
            if (Ls.eq.1) then
                pixo(i,j) = p_l(1)
                Lu = 1
                goto 9
                endif
            if (Ls.eq.2) then
                if (p_l(1).gt.p_l(2)) then
                    temp = p_l(1)
                    p_l(1) = p_l(2)
                    p_l(2) = temp
                    endif
                pixo(i,j) = (p_l(1)+p_l(2))/2.0
                Lu = 2
                if (p_l(1).lt.0.and.p_l(2).gt.0) then
                    pixo(i,j) = p_l(2)
                    Lu = 1
                    goto 9
                    endif
                if (p_l(1).gt.0.and.
     .              p_l(2)-p_l(1).gt.5*sqrt(abs(p_l(1))+25)) then
                    pixo(i,j) = p_l(1)
                    Lu = 1
                    goto 9
                    endif
                endif
    9       continue
            if (ig.eq.5001.and.Ls.gt.0.and.jg/100*100.eq.jg) then
                write( *,114) V1I2,
     .                        i, j, ig, jg, Ls, Lu,
     .                        pixo(i,j), (p_l(L),L=1,Ls)
                endif
  114       format(1x,'V1I2: ',i1,
     .             1x,2i5,5x,2i5,5x,i1,1x,i1,
     .                    5x,f8.2,5x,f8.2,1x,f8.2)
            enddo!i
            enddo!j
         if (V1I2.eq.1) FILENAME = 'dex_no_gaia_STEP10A_F606W.fits'
         if (V1I2.eq.2) FILENAME = 'dex_no_gaia_STEP10B_F814W.fits'
         call writfits_r4_WCSLTV(FILENAME,pixo,UDIM,VDIM,
     .                              CRPIX1,CRPIX2,CRVAL1,CRVAL2,
     .                              CD1_1, CD1_2, CD2_1, CD2_2,
     .                              LTV1, LTV2,
     .                              LTM1_1,LTM1_2,LTM2_1,LTM2_2)
         enddo
 
      deallocate(pixo)
 
 
      open(71,file='dex_no_gaia_STEP10C_GAIA.reg',status='unknown')
      write(71,'(''# Region file format: DS9 version 4.1'')')
      write(71,'(''global color=yellow dashlist=8 3 width=3'')')
      write(71,'(''fk5'')')
      do G = 1, Gs
         uu_g(G) = 0
         sz_g(G) = 0.25
         enddo
      do U = 1, Us
         if (uu_u(U).eq.1) then
             um_g(g_u(U)) = 1
             sz_g(g_u(U)) = 0.75
             endif
         enddo
      do G = 1, Gs
         if (um_g(G).eq.0) write(71,711) ra_g(G), de_g(G), sz_g(G)
         if (um_g(G).eq.1) write(71,712) ra_g(G), de_g(G), sz_g(G)
  711    format('circle(',f12.8,',',f12.8,',',f4.2,'") # color=red')
  712    format('circle(',f12.8,',',f12.8,',',f4.2,'") # color=green')
         enddo
      close(71)
 
      open(72,file='dex_no_gaia_STEP10D_HST.reg',status='unknown')
      write(72,'(''# Region file format: DS9 version 4.1'')')
      write(72,'(''global color=yellow width=1'')')
      write(72,'(''physical'')')
      do Q = 1, Qs
         siz = 0.20
         if (mv_q(Q).lt. -6.0) siz = 0.30
         if (mv_q(Q).lt. -7.0) siz = 0.40
         if (mv_q(Q).lt. -8.0) siz = 0.50
         if (mv_q(Q).lt. -9.0) siz = 0.60
         if (mv_q(Q).lt.-10.0) siz = 0.70
         if (mv_q(Q).lt.-11.0) siz = 0.80
         if (mv_q(Q).lt.-12.0) siz = 0.90
         if (mv_q(Q).lt.-13.0) siz = 1.00
         if (mv_q(Q).lt.-14.0) siz = 1.50
         if (mv_q(Q).lt.-15.0) siz = 2.00
         if (mv_q(Q).lt.-16.0) siz = 2.50
         if (mv_q(Q).lt.-17.0) siz = 3.00
         if (mv_q(Q).lt.-18.0) siz = 3.50
         if (mv_q(Q).lt.-19.0) siz = 4.00
         if (mv_q(Q).lt.-20.0) siz = 4.90
         write(72,171) uu_q(Q), vv_q(Q), siz*2
  171    format('circle(',f8.3,',',f8.3,',',f4.2,')')
         enddo
      close(72)
 
      STOP 'at end'
      end
 
 
 
c /grp/webpages/jayander/HST1PASS/LIB/GDCs/STDGDCs/ACSWFC/STDGDC_OFFICIAL_JFRAME_ACSWFC_F814W.fits
c /grp/webpages/jayander/HST1PASS/LIB/GDCs/STDGDCs/WFC3UV/STDGDC_WFC3UV_F814W.fits
 
 
 
c-----------------------------------------
c
c /user/jayander/FORTRAN/ROUTINES/STRING/str_contains.f
c
c this routine will return whether STRING1
c contains an incidence STRING2
c
       logical function str_contains(STRING1,L1,STRING2,L2)
       implicit none
 
       character*(*) STRING1
       integer L1
       character*(*) STRING2
       integer L2
 
       integer L1u
       integer L2u
 
       integer L1t
 
       L1u = L1
       L2u = L2
 
       if (L1u.eq.0) then
           L1u = 0
   11      continue
           if (STRING1(L1u+1:L1u+1).eq.' ') goto 12
           goto 11
   12      continue
           endif
 
       if (L2u.eq.0) then
           L2u = 0
   21      continue
           if (STRING1(L2u+1:L2u+1).eq.' ') goto 22
           goto 21
   22      continue
           endif
 
       str_contains = .true.
 
       do L1t = 1, L1u-L2u+1
          if (STRING1(L1t:L1t+L2u-1).eq.STRING2(1:L2u)) return
          enddo
 
       str_contains = .false.
       return
       end
 
 
 
 
 
 
 
c
c /user/jayander/FORTRAN/./ROUTINES/BARSIG/dbarsigs_2d.f
c
      subroutine dbarsigs_2d(xlist,ylist,NTOT,
     .                       xbar,xsig,
     .                       ybar,ysig,NUSE,SIGCLIP)
      implicit none
 
      integer NTOT
      real*8 xlist(NTOT), ylist(NTOT)
      real*8 xbar, ybar
      real*8 xsig, ysig
      integer NUSE
      real SIGCLIP
 
      integer n
      integer nsum, nsumo
      integer NIT
      real dsig
 
      real*8 xsumb, xsume
      real*8 ysumb, ysume
 
      nsum = 0
      xbar = 0.e0
      ybar = 0.e0
      xsig = 9e9
      ysig = 9e9
      do NIT = 1, 20
         dsig = sqrt(xsig**2+ysig**2)
         nsumo = nsum
         xsumb = 0.
         xsume = 0.
         ysumb = 0.
         ysume = 0.
         nsum = 0.
         do n = 1, NTOT
            if (sqrt((xlist(n)-xbar)**2+(ylist(n)-ybar)**2)
     .                  .le.dsig*SIGCLIP) then
               xsumb = xsumb + xlist(n)
               ysumb = ysumb + ylist(n)
               xsume = xsume + abs(xlist(n)-xbar)
               ysume = ysume + abs(ylist(n)-ybar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) then
             xbar = xsumb / nsum
             ybar = ysumb / nsum
             endif
         if (nsum.gt.1) then
             xsig = xsume/(nsum-1)
             ysig = ysume/(nsum-1)
             endif
        !print*,'---> NIT: ',NIT,xbar,ybar,xsig,ysig,dsig,NSUM
         if (nsum.lt.0.35*NTOT.and.NIT.ge.3) return
         if (nsum.eq.nsumo.and.NIT.ge.3) goto 1
         enddo
  1   continue
      NUSE = nsum
      if (nsum.le.1) then
          xsig = 0.999
          ysig = 0.999
          endif
 
      return
      end
 
 
c----------------------------------------------
c
c reports whether a point is inside a convex
c polygon
c
      logical function inside_poly(x,y,xl,yl,Ns)
      implicit none
 
      real*8 x,y
      integer Ns
      real*8 xl(Ns), yl(Ns)
 
      integer n
      real*8  cross1
      real*8  crossN
 
      inside_poly = .true.
      cross1 = (x-xl(Ns))*(yl(1)-yl(Ns))
     .       - (y-yl(Ns))*(xl(1)-xl(Ns))
      do n = 1, Ns-1
           crossN = (x-xl(N))*(yl(N+1)-yl(N))
     .            - (y-yl(N))*(xl(N+1)-xl(N))
           if (cross1*crossN.le.0) inside_poly = .false.
           enddo
      return
      end
 
 
