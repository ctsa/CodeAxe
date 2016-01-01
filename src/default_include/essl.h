/* Normal C Path   */
 
#ifndef __essl
 
  #define __essl 1
 
/**********************************************************************
*    LICENSED MATERIALS - PROPERTY OF IBM                             *
*    THIS MODULE IS "RESTRICTED MATERIALS OF IBM"                     *
*                                                                     *
*    IBM CONFIDENTIAL                                                 *
*                                                                     *
*    OCO SOURCE MATERIALS                                             *
*                                                                     *
*    5765-042                                                         *
*    5765-C42                                                         *
*                                                                     *
*    (C) COPYRIGHT IBM CORP. 1991, 2000.                              *
*                                                                     *
*    THE SOURCE CODE FOR THIS PROGRAM IS NOT PUBLISHED OR OTHERWISE   *
*    DIVESTED OF ITS TRADE SECRETS, IRRESPECTIVE OF WHAT HAS BEEN     *
*    DEPOSITED WITH THE U.S. COPYRIGHT OFFICE.                        *
*                                                                     *
*  Program name - <essl.h> header file                                *
*  Descriptive name - ESSL C,C++ language header file                 *
*                                                                     *
*  Function : This file must be included in any C file                *
*             containing ESSL calls in order for the ESSL calls       *
*             to work as documented in the ESSL Guide and             *
*  Change activity -                                                  *
*                    Added new routines for V2.2                      *
*                    Total = 425+16=441                (5/93)         *
*                    Added iessl()                      3/94          *
*                    Added extern for C++ path          6/94          *
*                    and fixed esvi4=i3 problem         6/94          *
*                    If CMPLX and DCMPLX not defined, defined  5/95   *
*                    Added CWLEV/ZWLEV and Check for aux1      5/96   *
*                    Rewrite for thread safe ESSL              1/97   *
*                    Change variable names that start with __  1/00   *
*                    to _ to fix a problem with C V5 Compiler         *
*                    Change the definition of the CMPLX C++ Class     *
*                    Added DBSTRF, DBSTRS, DBSSV and DGEQRF    5/00   *
**********************************************************************/
 
#ifndef  _CMPLX
#define  _CMPLX 1
#ifndef  _REIM
#define _REIM   1
#endif
typedef union { struct { float  _re, _im; } _data; double _align; }  cmplx;
#endif
 
#ifndef  _DCMPLX
#define  _DCMPLX 1
#ifndef  _REIM
#define _REIM   1
#endif
typedef union { struct { double _re, _im; } _data; double _align; } dcmplx;
#endif
 
#ifdef  _REIM
#define RE(x)   ((x)._data._re)
#define IM(x)   ((x)._data._im)
#endif
 
/*  Linear Algebra Subprograms  */
 
/*  Vector-Scalar Subprograms  */
 
#define isamax    esvisamax
#define idamax    esvidamax
#define icamax    esvicamax
#define izamax    esvizamax
 
#define isamin    esvisamin
#define idamin    esvidamin
 
#define ismax     esvismax
#define idmax     esvidmax
 
#define ismin     esvismin
#define idmin     esvidmin
 
#define sasum     esvsasum
#define dasum     esvdasum
#define scasum    esvscasum
#define dzasum    esvdzasum
 
#define saxpy     esvsaxpy
#define daxpy     esvdaxpy
#define caxpy     esvcaxpy
#define zaxpy     esvzaxpy
 
#define scopy     esvscopy
#define dcopy     esvdcopy
#define ccopy     esvccopy
#define zcopy     esvzcopy
 
#define sdot      esvsdot
#define ddot      esvddot
 
#define cdotu     esvcdotu
#define zdotu     esvzdotu
#define cdotc     esvcdotc
#define zdotc     esvzdotc
 
#define snaxpy    esvsnaxpy
#define dnaxpy    esvdnaxpy
#define sndot     esvsndot
#define dndot     esvdndot
 
#define snrm2     esvsnrm2
#define dnrm2     esvdnrm2
#define scnrm2    esvscnrm2
#define dznrm2    esvdznrm2
 
#define snorm2    esvsnorm2
#define dnorm2    esvdnorm2
#define cnorm2    esvcnorm2
#define znorm2    esvznorm2
 
#define srotg     esvsrotg
#define drotg     esvdrotg
#define crotg     esvcrotg
#define zrotg     esvzrotg
 
#define srot      esvsrot
#define drot      esvdrot
#define crot      esvcrot
#define zrot      esvzrot
 
#define csrot     esvcsrot
#define zdrot     esvzdrot
 
#define sscal     esvsscal
#define dscal     esvdscal
#define cscal     esvcscal
#define zscal     esvzscal
 
#define csscal    esvcsscal
#define zdscal    esvzdscal
 
#define sswap     esvsswap
#define dswap     esvdswap
#define cswap     esvcswap
#define zswap     esvzswap
 
#define syax      esvsyax
#define dyax      esvdyax
#define cyax      esvcyax
#define zyax      esvzyax
 
#define csyax     esvcsyax
#define zdyax     esvzdyax
 
#define szaxpy    esvszaxpy
#define dzaxpy    esvdzaxpy
#define czaxpy    esvczaxpy
#define zzaxpy    esvzzaxpy
 
#define svea      esvsvea
#define dvea      esvdvea
#define cvea      esvcvea
#define zvea      esvzvea
 
#define sves      esvsves
#define dves      esvdves
#define cves      esvcves
#define zves      esvzves
 
#define svem      esvsvem
#define dvem      esvdvem
#define cvem      esvcvem
#define zvem      esvzvem
 
 
/*  Sparse Vector-Scalar Subroutines  */
 
#define ssctr      esvssctr
#define dsctr      esvdsctr
#define csctr      esvcsctr
#define zsctr      esvzsctr
 
#define sgthr      esvsgthr
#define dgthr      esvdgthr
#define cgthr      esvcgthr
#define zgthr      esvzgthr
 
#define sgthrz     esvsgthrz
#define dgthrz     esvdgthrz
#define cgthrz     esvcgthrz
#define zgthrz     esvzgthrz
 
#define saxpyi     esvsaxpyi
#define daxpyi     esvdaxpyi
#define caxpyi     esvcaxpyi
#define zaxpyi     esvzaxpyi
 
#define sdoti      esvsdoti
#define ddoti      esvddoti
 
#define cdotui     esvcdotui
#define zdotui     esvzdotui
 
#define cdotci     esvcdotci
#define zdotci     esvzdotci
 
/*  Dense Matrix-Vector Subroutines  */
 
#define sgemv      esvsgemv
#define dgemv      esvdgemv
#define cgemv      esvcgemv
#define zgemv      esvzgemv
 
#define sgemx      esvsgemx
#define dgemx      esvdgemx
 
#define sgemtx     esvsgemtx
#define dgemtx     esvdgemtx
 
#define sger       esvsger1
#define dger       esvdger1
#define sger1      esvsger1
#define dger1      esvdger1
 
#define cgeru      esvcgeru
#define zgeru      esvzgeru
 
#define cgerc      esvcgerc
#define zgerc      esvzgerc
 
#define sslmx      esvsslmx
#define dslmx      esvdslmx
 
#define sslr1      esvsslr1
#define dslr1      esvdslr1
 
#define sslr2      esvsslr2
#define dslr2      esvdslr2
 
#define strmv      esvstrmv
#define dtrmv      esvdtrmv
#define ctrmv      esvctrmv
#define ztrmv      esvztrmv
 
#define sspmv      esvsspmv
#define dspmv      esvdspmv
#define chpmv      esvchpmv
#define zhpmv      esvzhpmv
 
#define ssymv      esvssymv
#define dsymv      esvdsymv
#define chemv      esvchemv
#define zhemv      esvzhemv
 
#define sspr       esvsspr
#define dspr       esvdspr
#define chpr       esvchpr
#define zhpr       esvzhpr
 
#define ssyr       esvssyr
#define dsyr       esvdsyr
#define cher       esvcher
#define zher       esvzher
 
#define sspr2      esvsspr2
#define dspr2      esvdspr2
#define chpr2      esvchpr2
#define zhpr2      esvzhpr2
 
#define ssyr2      esvssyr2
#define dsyr2      esvdsyr2
#define cher2      esvcher2
#define zher2      esvzher2
 
#define sgbmv      esvsgbmv
#define dgbmv      esvdgbmv
#define cgbmv      esvcgbmv
#define zgbmv      esvzgbmv
 
#define ssbmv      esvssbmv
#define dsbmv      esvdsbmv
#define chbmv      esvchbmv
#define zhbmv      esvzhbmv
 
#define stbmv      esvstbmv
#define dtbmv      esvdtbmv
#define ctbmv      esvctbmv
#define ztbmv      esvztbmv
 
#define stpmv      esvstpmv
#define dtpmv      esvdtpmv
#define ctpmv      esvctpmv
#define ztpmv      esvztpmv
 
/*  Sparse Matrix-Vector Linear Algebra Subroutines  */
 
 
#define dsmmx      esvdsmmx
 
#ifdef  __ESVERR
#define dsmtm      esvdsmtm_er
#else
#define dsmtm      esvdsmtm
#endif
 
#define dsdmx      esvdsdmx
 
/*  Matrix Operation Subroutines  */
 
#define sgeadd     esvsgeadd
#define dgeadd     esvdgeadd
#define cgeadd     esvcgeadd
#define zgeadd     esvzgeadd
 
#define sgesub     esvsgesub
#define dgesub     esvdgesub
#define cgesub     esvcgesub
#define zgesub     esvzgesub
 
#define sgemul     esvsgemul
#define dgemul     esvdgemul
#define cgemul     esvcgemul
#define zgemul     esvzgemul
 
#ifdef  __ESVERR
#define sgemms     esvsgemms_er
#define dgemms     esvdgemms_er
#define cgemms     esvcgemms_er
#define zgemms     esvzgemms_er
#else
#define sgemms     esvsgemms
#define dgemms     esvdgemms
#define cgemms     esvcgemms
#define zgemms     esvzgemms
#endif
 
#define sgemm      esvsgemm
#define dgemm      esvdgemm
#define cgemm      esvcgemm
#define zgemm      esvzgemm
 
#define ssyrk      esvssyrk
#define dsyrk      esvdsyrk
 
#define sgetmi     esvsgetmi
#define dgetmi     esvdgetmi
#define cgetmi     esvcgetmi
#define zgetmi     esvzgetmi
 
#define sgetmo     esvsgetmo
#define dgetmo     esvdgetmo
#define cgetmo     esvcgetmo
#define zgetmo     esvzgetmo
 
#define strmm      esvstrmm
#define dtrmm      esvdtrmm
#define ctrmm      esvctrmm
#define ztrmm      esvztrmm
 
#define ssymm      esvssymm
#define dsymm      esvdsymm
#define csymm      esvcsymm
#define zsymm      esvzsymm
 
#define ssyr2k     esvssyr2k
#define dsyr2k     esvdsyr2k
#define csyr2k     esvcsyr2k
#define zsyr2k     esvzsyr2k
 
#define csyrk      esvcsyrk
#define zsyrk      esvzsyrk
 
#define cherk      esvcherk
#define zherk      esvzherk
 
#define cher2k     esvcher2k
#define zher2k     esvzher2k
 
#define chemm      esvchemm
#define zhemm      esvzhemm
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
#define sgef       esvsgef
#define dgef       esvdgef
#define cgef       esvcgef
#define zgef       esvzgef
 
#define sges       esvsges
#define dges       esvdges
#define cges       esvcges
#define zges       esvzges
 
#define sgesm      esvsgesm
#define dgesm      esvdgesm
#define cgesm      esvcgesm
#define zgesm      esvzgesm
 
#ifdef  __ESVERR
#define sgefcd     esvsgefcd_er
#define dgefcd     esvdgefcd_er
#else
#define sgefcd     esvsgefcd
#define dgefcd     esvdgefcd
#endif
 
#define sppf       esvsppf
#define dppf       esvdppf
 
#define spps       esvspps
#define dpps       esvdpps
 
#ifdef  __ESVERR
#define sppfcd     esvsppfcd_er
#define dppfcd     esvdppfcd_er
 
#define sgeicd     esvsgeicd_er
#define dgeicd     esvdgeicd_er
 
#define sppicd     esvsppicd_er
#define dppicd     esvdppicd_er
 
#else
#define sppfcd     esvsppfcd
#define dppfcd     esvdppfcd
 
#define sgeicd     esvsgeicd
#define dgeicd     esvdgeicd
 
#define sppicd     esvsppicd
#define dppicd     esvdppicd
 
#endif
 
#define strsm      esvstrsm
#define dtrsm      esvdtrsm
#define ctrsm      esvctrsm
#define ztrsm      esvztrsm
 
#define strsv      esvstrsv
#define dtrsv      esvdtrsv
#define ctrsv      esvctrsv
#define ztrsv      esvztrsv
 
#define spof       esvspof
#define dpof       esvdpof
#define cpof       esvcpof
#define zpof       esvzpof
 
#define sposm      esvsposm
#define dposm      esvdposm
#define cposm      esvcposm
#define zposm      esvzposm
 
#ifdef  __ESVERR
#define spofcd     esvspofcd_er
#define dpofcd     esvdpofcd_er
 
#define spoicd     esvspoicd_er
#define dpoicd     esvdpoicd_er
 
#else
#define spofcd     esvspofcd
#define dpofcd     esvdpofcd
 
#define spoicd     esvspoicd
#define dpoicd     esvdpoicd
#endif
 
#define stri       esvstri
#define dtri       esvdtri
 
#define stpi       esvstpi
#define dtpi       esvdtpi
 
#define stpsv      esvstpsv
#define dtpsv      esvdtpsv
#define ctpsv      esvctpsv
#define ztpsv      esvztpsv
 
#define stbsv      esvstbsv
#define dtbsv      esvdtbsv
#define ctbsv      esvctbsv
#define ztbsv      esvztbsv
 
#define cgtnpf     esvcgtnpf
#define zgtnpf     esvzgtnpf
 
#define cgtnps     esvcgtnps
#define zgtnps     esvzgtnps
 
#define sgetrf     esvsgetrf
#define dgetrf     esvdgetrf
#define cgetrf     esvcgetrf
#define zgetrf     esvzgetrf
 
#define sgetrs     esvsgetrs
#define dgetrs     esvdgetrs
#define cgetrs     esvcgetrs
#define zgetrs     esvzgetrs
#define dbstrs     esvdbstrs
#ifdef  __ESVERR
#define dbstrf     esvdbstrf_er
#define dbssv      esvdbssv_er
#else
#define dbstrf     esvdbstrf
#define dbssv      esvdbssv
#endif
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
#define sgbf       esvsgbf
#define dgbf       esvdgbf
 
#define sgbs       esvsgbs
#define dgbs       esvdgbs
 
#define spbf       esvspbf
#define dpbf       esvdpbf
 
#define spbs       esvspbs
#define dpbs       esvdpbs
 
#define spbchf     esvspbchf
#define dpbchf     esvdpbchf
 
#define spbchs     esvspbchs
#define dpbchs     esvdpbchs
 
#define sgtf       esvsgtf
#define dgtf       esvdgtf
 
#define sgts       esvsgts
#define dgts       esvdgts
 
#define sgtnpf     esvsgtnpf
#define dgtnpf     esvdgtnpf
 
#define sgtnps     esvsgtnps
#define dgtnps     esvdgtnps
 
#define sgtnp      esvsgtnp
#define dgtnp      esvdgtnp
#define cgtnp      esvcgtnp
#define zgtnp      esvzgtnp
 
#define sptf       esvsptf
#define dptf       esvdptf
 
#define spts       esvspts
#define dpts       esvdpts
 
/* Sparse Linear Algebraic Equations Subroutines   */
 
#ifdef  __ESVERR
#define dgsf       esvdgsf_er
#define dgss       esvdgss_er
#define dgkfs      esvdgkfs_er
#define dskfs      esvdskfs_er
#define dsmcg      esvdsmcg_er
#define dsdcg      esvdsdcg_er
#define dsmgcg     esvdsmgcg_er
#define dsdgcg     esvdsdgcg_er
#else
#define dgsf       esvdgsf
#define dgss       esvdgss
#define dgkfs      esvdgkfs
#define dskfs      esvdskfs
#define dsmcg      esvdsmcg
#define dsdcg      esvdsdcg
#define dsmgcg     esvdsmgcg
#define dsdgcg     esvdsdgcg
#endif
 
#ifdef  __ESVERR
#define dsris      esvdsris_er
#else
#define dsris      esvdsris
#endif
 
/* Linear Least Square Subroutines  */
 
#ifdef  __ESVERR
#define sgesvf     esvsgesvf_er
#define dgesvf     esvdgesvf_er
#define sgells     esvsgells_er
#define dgells     esvdgells_er
#else
#define sgesvf     esvsgesvf
#define dgesvf     esvdgesvf
#define sgells     esvsgells
#define dgells     esvdgells
#endif
#define sgesvs     esvsgesvs
#define dgesvs     esvdgesvs
#define dgeqrf     esvdgeqrf
 
/* Eigensystem Analysis   Subroutines   */
 
#ifdef  __ESVERR
#define sgeev      esvsgeev_er
#define dgeev      esvdgeev_er
#define cgeev      esvcgeev_er
#define zgeev      esvzgeev_er
#define sslev  esvsspev_er
#define dslev  esvdspev_er
#define chlev  esvchpev_er
#define zhlev  esvzhpev_er
#define sspev      esvsspev_er
#define dspev      esvdspev_er
#define chpev      esvchpev_er
#define zhpev      esvzhpev_er
#define sspsv      esvsspsv_er
#define dspsv      esvdspsv_er
#define chpsv      esvchpsv_er
#define zhpsv      esvzhpsv_er
#define sgegv      esvsgegv_er
#define dgegv      esvdgegv_er
#define ssygv      esvssygv_er
#define dsygv      esvdsygv_er
#else
#define sgeev      esvsgeev
#define dgeev      esvdgeev
#define cgeev      esvcgeev
#define zgeev      esvzgeev
#define sslev  esvsspev
#define dslev  esvdspev
#define chlev  esvchpev
#define zhlev  esvzhpev
#define sspev      esvsspev
#define dspev      esvdspev
#define chpev      esvchpev
#define zhpev      esvzhpev
#define sspsv      esvsspsv
#define dspsv      esvdspsv
#define chpsv      esvchpsv
#define zhpsv      esvzhpsv
#define sgegv      esvsgegv
#define dgegv      esvdgegv
#define ssygv      esvssygv
#define dsygv      esvdsygv
#endif
 
/*  Fourier Transforms  Subroutines   */
 
#ifdef __ESVERR
#define scft       esvscft_er
#define dcft       esvdcft_er
#define srcft      esvsrcft_er
#define drcft      esvdrcft_er
#define scrft      esvscrft_er
#define dcrft      esvdcrft_er
#define scft2      esvscft2_er
#define srcft2     esvsrcft2_er
#define scrft2     esvscrft2_er
#define dcft2      esvdcft2_er
#define drcft2     esvdrcft2_er
#define dcrft2     esvdcrft2_er
#define scft3      esvscft3_er
#define srcft3     esvsrcft3_er
#define scrft3     esvscrft3_er
#define dcft3      esvdcft3_er
#define drcft3     esvdrcft3_er
#define dcrft3     esvdcrft3_er
#define scosft     esvscosft_er
#define scosf      esvscosf_er
#define dcosf      esvdcosf_er
#define ssinf      esvssinf_er
#define dsinf      esvdsinf_er
#else
#define scft       esvscft
#define dcft       esvdcft
#define srcft      esvsrcft
#define drcft      esvdrcft
#define scrft      esvscrft
#define dcrft      esvdcrft
#define scft2      esvscft2
#define srcft2     esvsrcft2
#define scrft2     esvscrft2
#define dcft2      esvdcft2
#define drcft2     esvdrcft2
#define dcrft2     esvdcrft2
#define scft3      esvscft3
#define srcft3     esvsrcft3
#define scrft3     esvscrft3
#define dcft3      esvdcft3
#define drcft3     esvdrcft3
#define dcrft3     esvdcrft3
#define scosft     esvscosft
#define scosf      esvscosf
#define dcosf      esvdcosf
#define ssinf      esvssinf
#define dsinf      esvdsinf
#endif
 
/*  Convulutions/Correlation Subroutines   */
 
#ifdef __ESVERR
#define scon       esvscon_er
#define scor       esvscor_er
#define sconf      esvsconf_er
#define scorf      esvscorf_er
#define sacor      esvsacor_er
#define sacorf     esvsacorf_er
#else
#define scon       esvscon
#define scor       esvscor
#define sconf      esvsconf
#define scorf      esvscorf
#define sacor      esvsacor
#define sacorf     esvsacorf
#endif
 
#define scond      esvscond
#define scord      esvscord
#define sdcon      esvsdcon
#define ddcon      esvddcon
 
#define sdcor      esvsdcor
#define ddcor      esvddcor
 
/*  Related Computations  Subroutines   */
 
#define spoly      esvspoly
#define dpoly      esvdpoly
 
#define sizc       esvsizc
#define dizc       esvdizc
 
#define strec      esvstrec
#define dtrec      esvdtrec
 
#define sqint      esvsqint
#define dqint      esvdqint
 
#ifdef __ESVERR
#define swlev      esvswlev_er
#define dwlev      esvdwlev_er
#define cwlev      esvcwlev_er
#define zwlev      esvzwlev_er
#else
#define swlev      esvswlev
#define dwlev      esvdwlev
#define cwlev      esvcwlev
#define zwlev      esvzwlev
#endif
 
/*  Sorting and Searching  Subroutines   */
 
#define isort      esvisort
#define ssort      esvssort
#define dsort      esvdsort
#define isortx     esvisortx
#define ssortx     esvssortx
#define dsortx     esvdsortx
 
#define ibsrch     esvibsrch
#define sbsrch     esvsbsrch
#define dbsrch     esvdbsrch
 
#define issrch     esvissrch
#define sssrch     esvsssrch
#define dssrch     esvdssrch
 
#define isorts     esvisorts
#define ssorts     esvssorts
#define dsorts     esvdsorts
 
/*  INTERPOLATION  Subroutines    */
 
#define spint      esvspint
#define dpint      esvdpint
 
#ifdef __ESVERR
#define stpint     esvstpint_er
#define dtpint     esvdtpint_er
#define scsin2     esvscsin2_er
#define dcsin2     esvdcsin2_er
#else
#define stpint     esvstpint
#define dtpint     esvdtpint
#define scsin2     esvscsin2
#define dcsin2     esvdcsin2
#endif
 
#define scsint     esvscsint
#define dcsint     esvdcsint
 
/*  NUMERICAL QUADRATURE  Subroutines    */
 
#define sptnq     esvsptnq
#define dptnq     esvdptnq
 
#define sglnq     esvsglnq
#define dglnq     esvdglnq
 
#define sglgq     esvsglgq
#define dglgq     esvdglgq
 
#define sgraq     esvsgraq
#define dgraq     esvdgraq
 
#define sghmq     esvsghmq
#define dghmq     esvdghmq
 
#define sglnq2    esvsglnq2
#define dglnq2    esvdglnq2
 
/*  Random Number Generator  Subroutines    */
 
#define surand    esvsurand
#define durand    esvdurand
 
#ifdef __ESVERR
#define snrand    esvsnrand_er
#define dnrand    esvdnrand_er
#else
#define snrand    esvsnrand
#define dnrand    esvdnrand
#endif
 
#define surxor    esvsurxor
#define durxor    esvdurxor
 
/*  Utility Subroutines     */
 
#define einfo     esveinfo
#define dsrsm     esvdsrsm
#define stride    esvstride
#define ivsset    esvivsset
#define ievops    esvievops
 
#ifdef __ESVERR
#define dgktrn    esvdgktrn_er
#define dsktrn    esvdsktrn_er
#else
#define dgktrn    esvdgktrn
#define dsktrn    esvdsktrn
#endif
 
/*  Parallel Processing Subroutines    */
 
#define dgemlp    esvdgemlp
 
#ifdef __ESVERR
#define dgefp     esvdgefp_er
#define dppfp     esvdppfp_er
#define dskfsp    esvdskfsp_er
#define dgkfsp    esvdgkfsp_er
#define scftp     esvscftp_er
#define scft2p    esvscft2p_er
#define scft3p    esvscft3p_er
#else
#define dgefp     esvdgefp
#define dppfp     esvdppfp
#define dskfsp    esvdskfsp
#define dgkfsp    esvdgkfsp
#define scftp     esvscftp
#define scft2p    esvscft2p
#define scft3p    esvscft3p
#endif
 
/* define fortran error routines name for consistent with main frame */
 
#define errset_   errset
#define errsav_   errsav
#define errstr_   errstr
 
 
/*  Linear Algebra Subprograms  */
 
/*  Vector-Scalar Subprograms  */
 
int esvisamax(int,  float *, int);
int esvidamax(int, double *, int);
int esvicamax(int,  cmplx *, int);
int esvizamax(int, dcmplx *, int);
 
int esvisamin(int,  float *, int);
int esvidamin(int, double *, int);
 
int esvismax(int,  float *, int);
int esvidmax(int, double *, int);
 
int esvismin(int,  float *, int);
int esvidmin(int, double *, int);
 
float  esvsasum(int,   float *, int);
double esvdasum(int,  double *, int);
float  esvscasum(int,  cmplx *, int);
double esvdzasum(int, dcmplx *, int);
 
void   esvsaxpy(int,  float,  float *, int,  float *, int);
void   esvdaxpy(int, double, double *, int, double *, int);
void   esvcaxpy(int,  cmplx,  cmplx *, int,  cmplx *, int);
void   esvzaxpy(int, dcmplx, dcmplx *, int, dcmplx *, int);
 
void   esvscopy(int,  float *, int,  float *, int);
void   esvdcopy(int, double *, int, double *, int);
void   esvccopy(int,  cmplx *, int,  cmplx *, int);
void   esvzcopy(int, dcmplx *, int, dcmplx *, int);
 
float  esvsdot(int,  float *, int,  float *, int);
double esvddot(int, double *, int, double *, int);
 
cmplx  esvcdotu(int,  cmplx *, int,  cmplx *, int);
dcmplx esvzdotu(int, dcmplx *, int, dcmplx *, int);
cmplx  esvcdotc(int, cmplx  *, int, cmplx  *, int);
dcmplx esvzdotc(int, dcmplx *, int, dcmplx *, int);
 
 
void   esvsnaxpy(int,int, float*, int,void *, int, int, void *, int, int);
void   esvdnaxpy(int,int,double*, int,void *, int, int, void *, int, int);
 
void   esvsndot(int,int, float *,int,int,void *,int,int,void *,int,int);
void   esvdndot(int,int,double *,int,int,void *,int,int,void *,int,int);
 
float  esvsnrm2(int,   float *, int);
double esvdnrm2(int,  double *, int);
float  esvscnrm2(int,  cmplx *, int);
double esvdznrm2(int, dcmplx *, int);
 
float  esvsnorm2(int,  float *, int);
double esvdnorm2(int, double *, int);
float  esvcnorm2(int,  cmplx *, int);
double esvznorm2(int, dcmplx *, int);
 
void   esvsrotg( float *,  float *,  float *,  float *);
void   esvdrotg(double *, double *, double *, double *);
void   esvcrotg( cmplx *,  cmplx *,  float *,  cmplx *);
void   esvzrotg(dcmplx *, dcmplx *, double *, dcmplx *);
 
void    esvsrot(int,  float *, int,  float *, int,  float,  float);
void    esvdrot(int, double *, int, double *, int, double, double);
void    esvcrot(int,  cmplx *, int,  cmplx *, int,  float,  cmplx);
void    esvzrot(int, dcmplx *, int, dcmplx *, int, double, dcmplx);
void   esvcsrot(int,  cmplx *, int,  cmplx *, int,  float,  float);
void   esvzdrot(int, dcmplx *, int, dcmplx *, int, double, double);
 
void    esvsscal(int,  float,  float *, int);
void    esvdscal(int, double, double *, int);
void    esvcscal(int,  cmplx,  cmplx *, int);
void    esvzscal(int, dcmplx, dcmplx *, int);
void   esvcsscal(int,  float,  cmplx *, int);
void   esvzdscal(int, double, dcmplx *, int);
 
void   esvsswap(int,  float *, int,  float *, int);
void   esvdswap(int, double *, int, double *, int);
void   esvcswap(int,  cmplx *, int,  cmplx *, int);
void   esvzswap(int, dcmplx *, int, dcmplx *, int);
 
void    esvsyax(int,  float,  float *, int,  float *, int);
void    esvdyax(int, double, double *, int, double *, int);
void    esvcyax(int,  cmplx,  cmplx *, int,  cmplx *, int);
void    esvzyax(int, dcmplx, dcmplx *, int, dcmplx *, int);
void   esvcsyax(int,  float,  cmplx *, int,  cmplx *, int);
void   esvzdyax(int, double, dcmplx *, int, dcmplx *, int);
 
void   esvszaxpy(int,  float,  float *, int,  float *, int,  float *, int);
void   esvdzaxpy(int, double, double *, int, double *, int, double *, int);
void   esvczaxpy(int,  cmplx,  cmplx *, int,  cmplx *, int,  cmplx *, int);
void   esvzzaxpy(int, dcmplx, dcmplx *, int, dcmplx *, int, dcmplx *, int);
 
void esvsvea(int,  float *, int,  float *, int,  float *, int);
void esvdvea(int, double *, int, double *, int, double *, int);
void esvcvea(int,  cmplx *, int,  cmplx *, int,  cmplx *, int);
void esvzvea(int, dcmplx *, int, dcmplx *, int, dcmplx *, int);
 
void esvsves(int,  float *, int,  float *, int,  float *, int);
void esvdves(int, double *, int, double *, int, double *, int);
void esvcves(int,  cmplx *, int,  cmplx *, int,  cmplx *, int);
void esvzves(int, dcmplx *, int, dcmplx *, int, dcmplx *, int);
 
void esvsvem(int,  float *, int,  float *, int,  float *, int);
void esvdvem(int, double *, int, double *, int, double *, int);
void esvcvem(int,  cmplx *, int,  cmplx *, int,  cmplx *, int);
void esvzvem(int, dcmplx *, int, dcmplx *, int, dcmplx *, int);
 
/*  Sparse Vector-Scalar Subroutines  */
 
void   esvssctr(int,  float *, int *,  float *);
void   esvdsctr(int, double *, int *, double *);
void   esvcsctr(int,  cmplx *, int *,  cmplx *);
void   esvzsctr(int, dcmplx *, int *, dcmplx *);
 
void   esvsgthr(int,  float *,  float *, int *);
void   esvdgthr(int, double *, double *, int *);
void   esvcgthr(int,  cmplx *,  cmplx *, int *);
void   esvzgthr(int, dcmplx *, dcmplx *, int *);
 
void   esvsgthrz(int,  float *,  float *, int *);
void   esvdgthrz(int, double *, double *, int *);
void   esvcgthrz(int,  cmplx *,  cmplx *, int *);
void   esvzgthrz(int, dcmplx *, dcmplx *, int *);
 
void   esvsaxpyi(int,  float,  float *, int *,  float *);
void   esvdaxpyi(int, double, double *, int *, double *);
void   esvcaxpyi(int,  cmplx,  cmplx *, int *,  cmplx *);
void   esvzaxpyi(int, dcmplx, dcmplx *, int *, dcmplx *);
 
float  esvsdoti(int,  float *, int *,  float *);
double esvddoti(int, double *, int *, double *);
 
cmplx  esvcdotui(int,  cmplx *, int *,  cmplx *);
dcmplx esvzdotui(int, dcmplx *, int *, dcmplx *);
 
cmplx  esvcdotci(int,  cmplx *, int *,  cmplx *);
dcmplx esvzdotci(int, dcmplx *, int *, dcmplx *);
 
/*  Dense Matrix-Vector Subroutines  */
 
void   esvsgemv(char *,int, int,  float, void *, int,  float *, int,  float,
                 float *, int);
void   esvdgemv(char *,int, int, double, void *, int, double *, int, double,
                double *, int);
void   esvcgemv(char *,int, int,  cmplx, void *, int,  cmplx *, int,  cmplx,
                 cmplx *, int);
void   esvzgemv(char *,int, int, dcmplx, void *, int, dcmplx *, int, dcmplx,
                dcmplx *, int);
 
void   esvsgemx(int,int, float,void *,int,  float *,int,  float *,int);
void   esvdgemx(int,int,double,void *,int, double *,int, double *,int);
 
void   esvsgemtx(int,int, float,void *,int, float *,int,  float *,int);
void   esvdgemtx(int,int,double,void *,int,double *,int, double *,int);
 
void   esvsger1(int, int,  float,  float *, int,  float *, int, void *, int);
void   esvdger1(int, int, double, double *, int, double *, int, void *, int);
void   esvcgeru(int, int,  cmplx,  cmplx *, int,  cmplx *, int, void *, int);
void   esvzgeru(int, int, dcmplx, dcmplx *, int, dcmplx *, int, void *, int);
void   esvcgerc(int, int,  cmplx,  cmplx *, int,  cmplx *, int, void *, int);
void   esvzgerc(int, int, dcmplx, dcmplx *, int, dcmplx *, int, void *, int);
 
void   esvsslmx(int,  float,  float *,  float *, int,  float *, int);
void   esvdslmx(int, double, double *, double *, int, double *, int);
 
void   esvsslr1(int,  float,  float *, int,  float *);
void   esvdslr1(int, double, double *, int, double *);
 
void   esvsslr2(int,  float,  float *, int,  float *, int,  float *);
void   esvdslr2(int, double, double *, int, double *, int, double *);
 
void esvstrmv(char *, char *, char *, int, void *, int,  float *, int);
void esvdtrmv(char *, char *, char *, int, void *, int, double *, int);
void esvctrmv(char *, char *, char *, int, void *, int,  cmplx *, int);
void esvztrmv(char *, char *, char *, int, void *, int, dcmplx *, int);
 
void esvsspmv(char *, int,  float,  float *,  float *, int,  float,  float *,
              int);
void esvdspmv(char *, int, double, double *, double *, int, double, double *,
              int);
void esvchpmv(char *, int,  cmplx,  cmplx *,  cmplx *, int,  cmplx,  cmplx *,
              int);
void esvzhpmv(char *, int, dcmplx, dcmplx *, dcmplx *, int, dcmplx, dcmplx *,
              int);
 
void esvssymv(char *, int,  float, void *, int,  float *, int,  float,
               float *, int);
void esvdsymv(char *, int, double, void *, int, double *, int, double,
              double *, int);
void esvchemv(char *, int,  cmplx, void *, int,  cmplx *, int,  cmplx,
               cmplx *, int);
void esvzhemv(char *, int, dcmplx, void *, int, dcmplx *, int, dcmplx,
              dcmplx *, int);
 
void esvsspr(char *, int,  float,  float *, int,   float *);
void esvdspr(char *, int, double, double *, int,  double *);
void esvchpr(char *, int,  float,  cmplx *, int,   cmplx *);
void esvzhpr(char *, int, double, dcmplx *, int,  dcmplx *);
 
void esvssyr(char *, int,  float,  float *, int,  void *, int);
void esvdsyr(char *, int, double, double *, int,  void *, int);
void esvcher(char *, int,  float,  cmplx *, int,  void *, int);
void esvzher(char *, int, double, dcmplx *, int,  void *, int);
 
void esvsspr2(char *, int,  float,  float *, int, float *, int,  float *);
void esvdspr2(char *, int, double, double *, int,double *, int, double *);
void esvchpr2(char *, int,  cmplx,  cmplx *, int, cmplx *, int,  cmplx *);
void esvzhpr2(char *, int, dcmplx, dcmplx *, int,dcmplx *, int, dcmplx *);
 
void esvssyr2(char *, int,  float,  float *, int,  float *, int, void *, int);
void esvdsyr2(char *, int, double, double *, int, double *, int, void *, int);
void esvcher2(char *, int,  cmplx,  cmplx *, int,  cmplx *, int, void *, int);
void esvzher2(char *, int, dcmplx, dcmplx *, int, dcmplx *, int, void *, int);
 
void esvsgbmv(char *, int, int, int, int,  float, void *, int,  float *, int,
                 float,  float *, int);
void esvdgbmv(char *, int, int, int, int, double, void *, int, double *, int,
                double, double *, int);
void esvcgbmv(char *, int, int, int, int,  cmplx, void *, int,  cmplx *, int,
                 cmplx,  cmplx *, int);
void esvzgbmv(char *, int, int, int, int, dcmplx, void *, int, dcmplx *, int,
                dcmplx, dcmplx *, int);
 
void esvssbmv(char *, int, int,  float, void *, int,  float *, int,  float,
               float *, int);
void esvdsbmv(char *, int, int, double, void *, int, double *, int, double,
              double *, int);
void esvchbmv(char *, int, int,  cmplx, void *, int,  cmplx *, int,  cmplx,
               cmplx *, int);
void esvzhbmv(char *, int, int, dcmplx, void *, int, dcmplx *, int, dcmplx,
              dcmplx *, int);
 
void esvstbmv(char *, char *, char *, int, int, void *, int, float *, int);
void esvdtbmv(char *, char *, char *, int, int, void *, int,double *, int);
void esvctbmv(char *, char *, char *, int, int, void *, int, cmplx *, int);
void esvztbmv(char *, char *, char *, int, int, void *, int,dcmplx *, int);
 
void esvstpmv(char *, char *, char *, int,  float *,  float *, int);
void esvdtpmv(char *, char *, char *, int, double *, double *, int);
void esvctpmv(char *, char *, char *, int,  cmplx *,  cmplx *, int);
void esvztpmv(char *, char *, char *, int, dcmplx *, dcmplx *, int);
 
/*  Sparse Matrix-Vector Subroutines  */
 
void  esvdsmmx(int, int, void *, void *, int, double *, double *);
 
#ifdef __ESVERR
 int  esvdsmtm_er(int, int, void *, void *, int, int *, int *, void *, void *,
                 int, float *, int *);
#else
 void esvdsmtm(int, int, void *, void *, int, int *, int *, void *, void *,
               int, float *, int);
#endif
 
void  esvdsdmx(int, int, int, void *, int, char *, int *, double *, double *);
 
/*  Matrix Operation Subroutines  */
 
void   esvsgeadd(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
void   esvdgeadd(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
void   esvcgeadd(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
void   esvzgeadd(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
 
void   esvsgesub(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
void   esvdgesub(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
void   esvcgesub(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
void   esvzgesub(void *, int, char *, void *, int, char *, void *, int, int,
                 int);
 
void   esvsgemul(void *, int, char *, void *, int, char *, void *, int, int,
                 int, int);
void   esvdgemul(void *, int, char *, void *, int, char *, void *, int, int,
                 int, int);
void   esvcgemul(void *, int, char *, void *, int, char *, void *, int, int,
                 int, int);
void   esvzgemul(void *, int, char *, void *, int, char *, void *, int, int,
                 int, int);
 
#ifdef __ESVERR
int esvsgemms_er(void *, int, char *, void *, int, char *, void *, int,
               int, int, int,  float *, int *);
int esvdgemms_er(void *, int, char *, void *, int, char *, void *, int,
               int, int, int, double *, int *);
int esvcgemms_er(void *, int, char *, void *, int, char *, void *, int,
               int, int, int,  float *, int *);
int esvzgemms_er(void *, int, char *, void *, int, char *, void *, int,
               int, int, int, double *, int *);
#else
void esvsgemms(void *, int, char *, void *, int, char *, void *, int,
               int, int, int,  float *, int);
void esvdgemms(void *, int, char *, void *, int, char *, void *, int,
               int, int, int, double *, int);
void esvcgemms(void *, int, char *, void *, int, char *, void *, int,
               int, int, int,  float *, int);
void esvzgemms(void *, int, char *, void *, int, char *, void *, int,
               int, int, int, double *, int);
#endif
 
void   esvsgemm(char *, char *, int, int, int,  float, void *, int, void *,
                int,  float, void *, int);
void   esvdgemm(char *, char *, int, int, int, double, void *, int, void *,
                int, double, void *, int);
void   esvcgemm(char *, char *, int, int, int,  cmplx, void *, int, void *,
                int,  cmplx, void *, int);
void   esvzgemm(char *, char *, int, int, int, dcmplx, void *, int, void *,
                int, dcmplx, void *, int);
 
void   esvssyrk(char *, char *, int, int,  float, void *, int,  float, void *,
                int);
void   esvdsyrk(char *, char *, int, int, double, void *, int, double, void *,
                int);
 
void   esvsgetmi(void *, int, int);
void   esvdgetmi(void *, int, int);
void   esvcgetmi(void *, int, int);
void   esvzgetmi(void *, int, int);
 
void   esvsgetmo(void *, int, int, int, void *, int);
void   esvdgetmo(void *, int, int, int, void *, int);
void   esvcgetmo(void *, int, int, int, void *, int);
void   esvzgetmo(void *, int, int, int, void *, int);
 
void esvstrmm(char *, char *, char *, char *, int, int, float, void*,
              int, void *, int);
void esvdtrmm(char *, char *, char *, char *, int, int, double, void*,
              int, void *, int);
void esvctrmm(char *, char *, char *, char *, int, int, cmplx, void*,
              int, void *, int);
void esvztrmm(char *, char *, char *, char *, int, int, dcmplx, void*,
              int, void *, int);
 
void esvssymm(char *, char *, int, int,  float, void *, int, void *, int,
               float, void *, int);
void esvdsymm(char *, char *, int, int, double, void *, int, void *, int,
              double, void *, int);
 
void esvssyr2k(char *, char *, int, int,  float, void *, int, void *, int,
                float, void *, int);
void esvdsyr2k(char *, char *, int, int, double, void *, int, void *, int,
               double, void *, int);
 
 
void esvcsyrk(char *, char *, int, int,  cmplx, void *, int,  cmplx, void *,
              int);
void esvzsyrk(char *, char *, int, int, dcmplx, void *, int, dcmplx, void *,
              int);
 
void esvcherk(char *, char *, int, int,  float, void *, int,  float, void *,
              int);
void esvzherk(char *, char *, int, int, double, void *, int, double, void *,
              int);
 
void esvcsyr2k(char *, char *, int, int,  cmplx, void *, int, void *, int,
                cmplx, void *, int);
void esvzsyr2k(char *, char *, int, int, dcmplx, void *, int, void *, int,
               dcmplx, void *, int);
void esvcher2k(char *, char *, int, int,  cmplx, void *, int, void *, int,
                float, void *, int);
void esvzher2k(char *, char *, int, int, dcmplx, void *, int, void *, int,
               double, void *, int);
 
 
void esvcsymm(char *, char *, int, int,  cmplx, void *, int, void *, int,
                cmplx, void *, int);
void esvzsymm(char *, char *, int, int, dcmplx, void *, int, void *, int,
               dcmplx, void *, int);
void esvchemm(char *, char *, int, int,  cmplx, void *, int, void *, int,
                cmplx, void *, int);
void esvzhemm(char *, char *, int, int, dcmplx, void *, int, void *, int,
               dcmplx, void *, int);
 
/*  Dense Linear Algebraic Equation Subroutines  */
 
int esvsgef(void *, int, int, int *);
int esvdgef(void *, int, int, int *);
int esvcgef(void *, int, int, int *);
int esvzgef(void *, int, int, int *);
 
void   esvsges(void *, int, int, int *,  float *, int);
void   esvdges(void *, int, int, int *, double *, int);
void   esvcges(void *, int, int, int *,  cmplx *, int);
void   esvzges(void *, int, int, int *, dcmplx *, int);
 
void esvsgesm(char *, void *, int, int, int *, void *, int, int);
void esvdgesm(char *, void *, int, int, int *, void *, int, int);
void esvcgesm(char *, void *, int, int, int *, void *, int, int);
void esvzgesm(char *, void *, int, int, int *, void *, int, int);
 
#ifdef __ESVERR
int esvsgefcd_er(void *, int, int, int *, int, float *, float *, float *,
                 int *);
int esvdgefcd_er(void *, int, int, int *, int, double *, double *, double *,
                 int *);
#else
void esvsgefcd(void *, int, int, int *, int, float *, float *, float *, int);
void esvdgefcd(void *, int, int, int *, int,double *,double *,double *, int);
#endif
 
int esvsppf( float *, int, int);
int esvdppf(double *, int, int);
 
void   esvspps( float *, int,  float *, int);
void   esvdpps(double *, int, double *, int);
 
#ifdef __ESVERR
int esvsppfcd_er( float *, int, int,  float *,  float *,  float *, int *);
int esvdppfcd_er(double *, int, int, double *, double *, double *, int *);
 
int esvsgeicd_er(void *, int, int, int,  float *,  float *,  float *, int *);
int esvdgeicd_er(void *, int, int, int, double *, double *, double *, int *);
 
int esvsppicd_er( float *, int, int,  float *,  float *,  float *, int *);
int esvdppicd_er(double *, int, int, double *, double *, double *, int *);
#else
void esvsppfcd( float *, int, int,  float *,  float *,  float *, int);
void esvdppfcd(double *, int, int, double *, double *, double *, int);
 
void esvsgeicd(void *, int, int, int,  float *,  float *,  float *, int);
void esvdgeicd(void *, int, int, int, double *, double *, double *, int);
 
void esvsppicd( float *, int, int,  float *,  float *,  float *, int);
void esvdppicd(double *, int, int, double *, double *, double *, int);
#endif
 
void  esvstrsm(char *, char *, char *, char *, int, int,  float, void *,
              int, void *, int);
void  esvdtrsm(char *, char *, char *, char *, int, int, double, void *,
              int, void *, int);
void  esvctrsm(char *, char *, char *, char *, int, int,  cmplx, void *,
               int, void *, int);
void  esvztrsm(char *, char *, char *, char *, int, int, dcmplx, void *,
               int, void *, int);
void  esvstrsv(char *, char *, char *, int, void *, int, void *, int);
void  esvdtrsv(char *, char *, char *, int, void *, int, void *, int);
void  esvctrsv(char *, char *, char *, int, void *, int, void *, int);
void  esvztrsv(char *, char *, char *, int, void *, int, void *, int);
 
 
void  esvspof(char *, void *, int, int);
void  esvdpof(char *, void *, int, int);
void  esvcpof(char *, void *, int, int);
void  esvzpof(char *, void *, int, int);
 
void  esvsposm(char *, void *, int, int, void *, int, int);
void  esvdposm(char *, void *, int, int, void *, int, int);
void  esvcposm(char *, void *, int, int, void *, int, int);
void  esvzposm(char *, void *, int, int, void *, int, int);
 
#ifdef __ESVERR
int esvspofcd_er(char*, void*, int, int, int,  float*,  float*,  float*, int *);
int esvdpofcd_er(char*, void*, int, int, int, double*, double*, double*, int *);
 
int esvspoicd_er(char*, void*, int, int, int,  float*,  float*,  float*, int *);
int esvdpoicd_er(char*, void*, int, int, int, double*, double*, double*, int *);
#else
void esvspofcd(char *, void *, int, int, int,  float*,  float*, float *, int);
void esvdpofcd(char *, void *, int, int, int, double*, double*, double*, int);
 
void esvspoicd(char *, void *, int, int, int,  float*,  float*, float *, int);
void esvdpoicd(char *, void *, int, int, int, double*, double*, double*, int);
#endif
 
int esvstri(char *, char *, void *, int, int);
int esvdtri(char *, char *, void *, int, int);
 
int esvstpi(char *, char *,  float *, int);
int esvdtpi(char *, char *, double *, int);
 
void esvstpsv(char *, char *, char *, int,  float *,  float *, int);
void esvdtpsv(char *, char *, char *, int, double *, double *, int);
void esvctpsv(char *, char *, char *, int,  cmplx *,  cmplx *, int);
void esvztpsv(char *, char *, char *, int, dcmplx *, dcmplx *, int);
 
void esvstbsv(char *, char *, char *, int, int, void *, int,  float *, int);
void esvdtbsv(char *, char *, char *, int, int, void *, int, double *, int);
void esvctbsv(char *, char *, char *, int, int, void *, int,  cmplx *, int);
void esvztbsv(char *, char *, char *, int, int, void *, int, dcmplx *, int);
 
void esvcgtnpf(int,  cmplx *,  cmplx *,  cmplx *, int);
void esvzgtnpf(int, dcmplx *, dcmplx *, dcmplx *, int);
 
void esvcgtnps(int,  cmplx *,  cmplx *,  cmplx *,  cmplx *);
void esvzgtnps(int, dcmplx *, dcmplx *, dcmplx *, dcmplx *);
 
void   esvsgetrf(int, int, void *, int,  int *, int *);
void   esvdgetrf(int, int, void *, int,  int *, int *);
void   esvcgetrf(int, int, void *, int,  int *, int *);
void   esvzgetrf(int, int, void *, int,  int *, int *);
 
void   esvsgetrs(char *,int,int,void *, int, int *, void *, int,int *);
void   esvdgetrs(char *,int,int,void *, int, int *, void *, int,int *);
void   esvcgetrs(char *,int,int,void *, int, int *, void *, int,int *);
void   esvzgetrs(char *,int,int,void *, int, int *, void *, int,int *);
void   esvdbstrs(char *, int, int, void *, int *, void *, int, int *);
#ifdef __ESVERR
int esvdbstrf_er(char *, int, void *, int *, int *);
int esvdbssv_er(char *, int, int, void *, int *, void *, int, int *);
#else
void esvdbstrf(char *, int, void *, int *, int *);
void esvdbssv(char *, int, int, void *, int *, void *, int, int *);
#endif
 
/*  Banded Linear Algebraic Equation Subroutines  */
 
int esvsgbf(void *, int, int, int, int, int *);
int esvdgbf(void *, int, int, int, int, int *);
 
void   esvsgbs(void *, int, int, int, int, int *,  float *);
void   esvdgbs(void *, int, int, int, int, int *, double *);
 
int esvspbf(void *, int, int, int);
int esvdpbf(void *, int, int, int);
 
void   esvspbs(void *, int, int, int,  float *);
void   esvdpbs(void *, int, int, int, double *);
 
int esvspbchf(void *, int, int, int);
int esvdpbchf(void *, int, int, int);
 
void   esvspbchs(void *, int, int, int,  float *);
void   esvdpbchs(void *, int, int, int, double *);
 
int esvsgtf(int,  float *,  float *,  float *,  float *, int *);
int esvdgtf(int, double *, double *, double *, double *, int *);
 
void   esvsgts(int,  float *,  float *,  float *,  float *, int *,  float *);
void   esvdgts(int, double *, double *, double *, double *, int *, double *);
 
void   esvsgtnpf(int,  float *,  float *,  float *, int);
void   esvdgtnpf(int, double *, double *, double *, int);
 
void   esvsgtnps(int,  float *,  float *,  float *,  float *);
void   esvdgtnps(int, double *, double *, double *, double *);
 
void   esvsgtnp(int,  float *,  float *,  float *,  float *);
void   esvdgtnp(int, double *, double *, double *, double *);
void   esvcgtnp(int,  cmplx *,  cmplx *,  cmplx *,  cmplx *);
void   esvzgtnp(int, dcmplx *, dcmplx *, dcmplx *, dcmplx *);
 
void   esvsptf(int,  float *,  float *, int);
void   esvdptf(int, double *, double *, int);
 
void   esvspts(int,  float *,  float *,  float *);
void   esvdpts(int, double *, double *, double *);
 
/*  Sparse Linear Algebraic Equation Subroutines  */
 
#ifdef __ESVERR
int esvdgsf_er(int, int, int, double *, int *, int *, int, int *, double *,
               double *, double *, int *);
int esvdgss_er(int, int, double *, int *, int *, int, double*, double*, int *);
 
int esvdgkfs_er(int, double *, int, int *, double *, int, int *, int *,
                double *, double *, int *, void *, int, int);
int esvdskfs_er(int, double *, int, int *, int *, double *, double *,
                int *, void *, int, int);
 
int esvdsmcg_er(int, int, void *, void *, int, double *, double *, int *,
                double *, double *, int *, double *, int *);
int esvdsdcg_er(int, int, int, void *, int, int *, double *, double *,
                int *, double *, double *, int *, double *, int *);
 
int esvdsmgcg_er(int, int, void *, void *, int, double *, double *, int *,
                 double *, double *, int *, double *, int *);
int esvdsdgcg_er(int, int, void *, int, int *, double *, double *, int *,
                 double *, double *, int *, double *, int *);
#else
void esvdgsf(int, int, int, double *, int *, int *, int, int *, double *,
             double *, double *, int);
void esvdgss(int, int, double *, int *, int *, int, double *, double *, int);
 
void esvdgkfs(int, double *, int, int *, double *, int, int *, int *,
              double *, double *, int, void *, int, int);
void esvdskfs(int, double *, int, int *, int *, double *, double *,
              int, void *, int, int);
 
void esvdsmcg(int, int, void *, void *, int, double *, double *, int *,
              double *, double *, int, double *, int);
void esvdsdcg(int, int, int, void *, int, int *, double *, double *,
              int *, double *, double *, int, double *, int);
 
void esvdsmgcg(int, int, void *, void *, int, double *, double *, int *,
               double *, double *, int, double *, int);
void esvdsdgcg(int, int, void *, int, int *, double *, double *, int *,
               double *, double *, int, double *, int);
#endif
 
#ifdef __ESVERR
int esvdsris_er(char *, char *, int, double *, int *, int *, double *,
                double *, int *, double *, double *, int *, double *, int *);
#else
void esvdsris(char *, char *, int, double *, int *, int *, double *,
              double *, int *, double *, double *, int, double *, int);
#endif
 
/*  Linear Least Squares Subroutines  */
 
#ifdef __ESVERR
int esvsgesvf_er(int, void *, int, void *, int, int,  float *, int,
                 int, float *, int *);
int esvdgesvf_er(int, void *, int, void *, int, int, double *, int,
                 int, double *, int *);
 
int esvsgells_er(int, void *, int, void *, int, void *, int,  float *,
                  float, int, int, int, int *,  float *, int *);
int esvdgells_er(int, void *, int, void *, int, void *, int, double *,
                 double, int, int, int, int *, double *, int *);
#else
void esvsgesvf(int, void *, int, void *, int, int,  float *, int,
               int,  float *, int);
void esvdgesvf(int, void *, int, void *, int, int, double *, int,
               int, double *, int);
 
void esvsgells(int, void *, int, void *, int, void *, int,  float *,
                float, int, int, int, int *,  float *, int);
void esvdgells(int, void *, int, void *, int, void *, int, double *,
               double, int, int, int, int *, double *, int);
#endif
 
void   esvsgesvs(void *, int, void *, int, int,  float *, void *, int,
                 int, int,  float);
void   esvdgesvs(void *, int, void *, int, int, double *, void *, int,
                 int, int, double);
 
void   esvdgeqrf(int, int, void *, int,  double *,double *, int, int *);
 
/*  Eigensystem Analysis Subroutines  */
 
#ifdef __ESVERR
int esvsgeev_er(int, void *, int,  cmplx *, void *, int, int *, int,
                 float *, int *);
int esvdgeev_er(int, void *, int, dcmplx *, void *, int, int *, int,
                double *, int *);
int esvcgeev_er(int, void *, int,  cmplx *, void *, int, int *, int,
                 float *, int *);
int esvzgeev_er(int, void *, int, dcmplx *, void *, int, int *, int,
                double *, int *);
#else
void esvsgeev(int, void *, int,  cmplx *, void *, int, int *, int,
               float *, int);
void esvdgeev(int, void *, int, dcmplx *, void *, int, int *, int,
              double *, int);
void esvcgeev(int, void *, int,  cmplx *, void *, int, int *, int,
               float *, int);
void esvzgeev(int, void *, int, dcmplx *, void *, int, int *, int,
              double *, int);
#endif
 
#ifdef __ESVERR
int esvsspev_er(int,  float *,  float *, void *, int, int,  float *, int *);
int esvdspev_er(int, double *, double *, void *, int, int, double *, int *);
int esvchpev_er(int,  cmplx *,  float *, void *, int, int,  float *, int *);
int esvzhpev_er(int, dcmplx *, double *, void *, int, int, double *, int *);
#else
void esvsspev(int,  float *,  float *, void *, int, int,  float *, int);
void esvdspev(int, double *, double *, void *, int, int, double *, int);
void esvchpev(int,  cmplx *,  float *, void *, int, int,  float *, int);
void esvzhpev(int, dcmplx *, double *, void *, int, int, double *, int);
#endif
 
#ifdef __ESVERR
int esvsspsv_er(int,  float *,  float *, void *, int, int, int,  float*, int*);
int esvdspsv_er(int, double *, double *, void *, int, int, int, double*, int*);
int esvchpsv_er(int,  cmplx *,  float *, void *, int, int, int,  float*, int*);
int esvzhpsv_er(int, dcmplx *, double *, void *, int, int, int, double*, int*);
#else
void esvsspsv(int,  float *,  float *, void *, int, int, int,  float*, int);
void esvdspsv(int, double *, double *, void *, int, int, int, double*, int);
void esvchpsv(int,  cmplx *,  float *, void *, int, int, int,  float*, int);
void esvzhpsv(int, dcmplx *, double *, void *, int, int, int, double*, int);
#endif
 
#ifdef __ESVERR
int esvsgegv_er(int, void *, int, void *, int,  cmplx *,  float *, void *,
                int, int,  float *, int *);
int esvdgegv_er(int, void *, int, void *, int, dcmplx *, double *, void *,
                int, int, double *, int *);
 
int esvssygv_er(int, void *, int, void *, int,  float *, void *, int,
                int, float *, int *);
int esvdsygv_er(int, void *, int, void *, int, double *, void *, int,
                int, double *, int *);
#else
void esvsgegv(int, void *, int, void *, int,  cmplx *,  float *, void *,
              int, int,  float *, int);
void esvdgegv(int, void *, int, void *, int, dcmplx *, double *, void *,
              int, int, double *, int);
 
void esvssygv(int, void *, int, void *, int,  float *, void *, int,
              int, float *, int);
void esvdsygv(int, void *, int, void *, int, double *, void *, int,
              int, double *, int);
#endif
 
/*  Fourier Transform Subroutines  */
 
#ifdef __ESVERR
int esvscft_er(int, void *, int, int, void *, int, int, int *, int,
               int, float, double *, int *, double *, int *);
int esvdcft_er(int, void *, int, int, void *, int, int, int *, int,
               int, double, double *, int *, double *, int *);
 
int esvsrcft_er(int, void *, int, void *, int, int *, int, int, float,
                double *, int *, double *, int *, double *, int *);
int esvdrcft_er(int, void *, int, void *, int, int *, int, int, double,
                double *, int *, double *, int *);
 
int esvscrft_er(int, void *, int, void *, int, int *, int, int, float,
                double *, int *, double *, int *, double *, int *);
 
int esvdcrft_er(int, void *, int, void *, int, int *, int, int, double,
                double *, int *, double *, int *);
#else
void esvscft(int, void *, int, int, void *, int, int, int, int,
             int, float, double *, int, double *, int);
void esvdcft(int, void *, int, int, void *, int, int, int, int,
             int, double, double *, int, double *, int);
 
void esvsrcft(int, void *, int, void *, int, int, int, int, float,
              double *, int, double *, int, double *, int);
void esvdrcft(int, void *, int, void *, int, int, int, int, double,
              double *, int, double *, int);
 
void esvscrft(int, void *, int, void *, int, int, int, int, float,
              double *, int, double *, int, double *, int);
 
void esvdcrft(int, void *, int, void *, int, int, int, int, double,
              double *, int, double *, int);
#endif
 
 
 
#ifdef __ESVERR
int esvdcft2_er(int, void *, int, int, void *, int, int, int *,
                int *, int, double, double *, int *, double *, int *);
 
int esvdrcft2_er(int, void *, int, void *, int, int *, int *, int,
                 double, double *, int *, double *, int *);
 
int esvdcrft2_er(int, void *, int, void *, int, int *, int *, int,
                 double, double *, int *, double *, int *);
 
int esvdcft3_er(void *, int, int, void *, int, int, int *, int *,
                int *, int, double, double *, int *);
 
int esvdrcft3_er(void *, int, int, void *, int, int, int *, int *,
                 int *, int, double, double *, int *);
 
int esvdcrft3_er(void *, int, int, void *, int, int, int *, int *,
                 int *, int, double, double *, int *);
 
int esvscosft_er(int, void *, int, int, void *, int, int, int *, int,
                 float, double *, int *, double *, int *);
 
int esvscosf_er(int, void *, int, int, void *, int, int, int *, int,
                 float, double *, int *, double *, int *);
 
int esvdcosf_er(int, void *, int, int, void *, int, int, int *, int,
                double, double *, int *, double *, int *);
 
int esvssinf_er(int, void *, int, int, void *, int, int, int *, int,
                float, double *, int *, double *, int *);
 
int esvdsinf_er(int, void *, int, int, void *, int, int, int *, int,
                double, double *, int *, double *, int *);
 
#else
void esvdcft2(int, void *, int, int, void *, int, int, int,
              int, int, double, double *, int, double *, int);
 
void esvdrcft2(int, void *, int, void *, int, int, int, int,
               double, double *, int, double *, int);
 
void esvdcrft2(int, void *, int, void *, int, int, int, int,
               double, double *, int, double *, int);
 
void esvdcft3(void *, int, int, void *, int, int, int, int,
              int, int, double, double *, int);
 
void esvdrcft3(void *, int, int, void *, int, int, int, int,
               int, int, double, double *, int);
 
void esvdcrft3(void *, int, int, void *, int, int, int, int,
               int, int, double, double *, int);
 
void esvscosft(int, void *, int, int, void *, int, int, int, int,
               float, double *, int, double *, int);
 
 
void esvscosf(int, void *, int, int, void *, int, int, int, int,
              float, double *, int, double *, int);
 
void esvdcosf(int, void *, int, int, void *, int, int, int, int,
              double, double *, int, double *, int);
 
void esvssinf(int, void *, int, int, void *, int, int, int, int,
              float, double *, int, double *, int);
 
void esvdsinf(int, void *, int, int, void *, int, int, int, int,
              double, double *, int, double *, int);
 
#endif
 
#ifdef __ESVERR
int esvscft2_er(int, void *, int, int, void *, int, int, int *, int *,
                int, float, double *, int *, double *, int *);
 
int esvsrcft2_er(int,void *, int, void *, int, int *, int *, int, float,
                 double *, int *, double *, int *, double *, int *);
 
int esvscrft2_er(int,void *, int, void *, int, int *, int *, int, float,
                 double *, int *, double *, int *, double *, int *);
 
int esvscft3_er(void *, int, int, void *, int, int, int *, int *, int *,
                int, float, double *, int *);
 
int esvsrcft3_er(void *, int, int, void *, int, int, int *, int *, int *,
                 int, float, double *, int *);
 
int esvscrft3_er(void *, int, int, void *, int, int, int *, int *, int *,
                 int, float, double *, int *);
#else
void esvscft2(int, void *, int, int, void *, int, int, int, int,
              int, float, double *, int, double *, int);
 
void esvsrcft2(int,void *, int, void *, int, int, int, int, float,
               double *, int, double *, int, double *, int);
 
void esvscrft2(int,void *, int, void *, int, int, int, int, float,
               double *, int, double *, int, double *, int);
 
void esvscft3(void *, int, int, void *, int, int, int, int, int,
              int, float, double *, int);
 
void esvsrcft3(void *, int, int, void *, int, int, int, int, int,
               int, float, double *, int);
 
void esvscrft3(void *, int, int, void *, int, int, int, int, int,
               int, float, double *, int);
#endif
 
/*  Convolutions/Correlations Subroutines  */
 
#ifdef __ESVERR
int esvscon_er(int, float *, int, void *, int, int, void *, int, int,
               int, int, int, int, int, double *, int *, double *, int *);
int esvscor_er(int, float *, int, void *, int, int, void *, int, int,
               int, int, int, int, int, double *, int *, double *, int *);
#else
void esvscon(int, float *, int, void *, int, int, void *, int, int,
             int, int, int, int, int, double *, int, double *, int);
void esvscor(int, float *, int, void *, int, int, void *, int, int,
             int, int, int, int, int, double *, int, double *, int);
#endif
 
void   esvscond(float *, int, void *, int, void *, int, int, int, int, int);
void   esvscord(float *, int, void *, int, void *, int, int, int, int, int);
 
#ifdef __ESVERR
int esvsconf_er(int, float *, int, void *, int, int, void *, int,
                int, int, int, int, int, int, double *, int *, double *, int *);
int esvscorf_er(int, float *, int, void *, int, int, void *, int,
                int, int, int, int, int, int, double *, int *, double *, int *);
#else
void esvsconf(int, float *, int, void *, int, int, void *, int,
              int, int, int, int, int, int, double *, int, double *, int);
void esvscorf(int, float *, int, void *, int, int, void *, int,
              int, int, int, int, int, int, double *, int, double *, int);
#endif
 
void   esvsdcon( float *, int, void *, int, void *, int, int, int, int, int,
                 int);
void   esvddcon(double *, int, void *, int, void *, int, int, int, int, int,
                 int);
void   esvsdcor( float *, int, void *, int, void *, int, int, int, int, int,
                 int);
void   esvddcor(double *, int, void *, int, void *, int, int, int, int, int,
                 int);
 
#ifdef __ESVERR
int esvsacor_er(int, void *, int, int, void *, int, int, int,
                int, int, double *, int *, double *, int *);
 
int esvsacorf_er(int, void *, int, int, void *, int, int, int,
                 int, int, double *, int *, double *, int *);
#else
void esvsacor(int, void *, int, int, void *, int, int, int,
              int, int, double *, int, double *, int);
 
void esvsacorf(int, void *, int, int, void *, int, int, int,
               int, int, double *, int, double *, int);
#endif
 
/*  Related Computations Subroutines  */
 
void   esvspoly( float *, int, int,  float *, int,  float *, int, int);
void   esvdpoly(double *, int, int, double *, int, double *, int, int);
 
void   esvsizc( float *, int, int, int, int *);
void   esvdizc(double *, int, int, int, int *);
 
void   esvstrec( float,  float *, int,  float *, int,  float *, int, int, int);
void   esvdtrec(double, double *, int, double *, int, double *, int, int, int);
 
int esvsqint( float,  float,  float,  float *, int, int,  float *, int,
              float *, int, int);
int esvdqint(double, double, double, double *, int, int, double *, int,
             double *, int, int);
 
#ifdef __ESVERR
int esvswlev_er( float *, int,  float *, int,  float *, int, int, double *,
                 int *);
int esvdwlev_er(double *, int, double *, int, double *, int, int, double *,
                 int *);
int esvcwlev_er( cmplx *, int,  cmplx *, int,  cmplx *, int, int, dcmplx *,
                 int *);
int esvzwlev_er(dcmplx *, int, dcmplx *, int, dcmplx *, int, int, dcmplx *,
                 int *);
#else
void esvswlev( float *, int,  float *, int,  float *, int, int, double *, int);
void esvdwlev(double *, int, double *, int, double *, int, int, double *, int);
void esvcwlev( cmplx *, int,  cmplx *, int,  cmplx *, int, int, dcmplx *, int);
void esvzwlev(dcmplx *, int, dcmplx *, int, dcmplx *, int, int, dcmplx *, int);
#endif
 
/*  Sorting and Searching Subroutines  */
 
void   esvisort(   int *, int, int);
void   esvssort( float *, int, int);
void   esvdsort(double *, int, int);
 
void   esvisortx(   int *, int, int, int *);
void   esvssortx( float *, int, int, int *);
void   esvdsortx(double *, int, int, int *);
 
void   esvibsrch(   int *, int, int,    int *, int, int, int *, int *, int);
void   esvsbsrch( float *, int, int,  float *, int, int, int *, int *, int);
void   esvdbsrch(double *, int, int, double *, int, int, int *, int *, int);
 
void   esvissrch(   int *, int, int,    int *, int, int, int, int *);
void   esvsssrch( float *, int, int,  float *, int, int, int, int *);
void   esvdssrch(double *, int, int, double *, int, int, int, int *);
 
void   esvisorts(   int *, int, int, int *,    int *, int);
void   esvssorts( float *, int, int, int *,  float *, int);
void   esvdsorts(double *, int, int, int *, double *, int);
 
/*  Interpolation Subroutines  */
 
void   esvspint( float *,  float *, int,  float *, int *,  float *,  float *,
                 int);
void   esvdpint(double *, double *, int, double *, int *, double *, double *,
                 int);
 
#ifdef __ESVERR
int esvstpint_er( float *,  float *, int, int,  float *,  float *, int,
                  float *, int *);
 
int esvdtpint_er(double *, double *, int, int, double *, double *, int,
                 double *, int *);
#else
void esvstpint( float *,  float *, int, int,  float *,  float *, int,
                float *, int);
 
void esvdtpint(double *, double *, int, int, double *, double *, int,
               double *, int);
#endif
 
void   esvscsint( float *,  float *, void *, int, int *,  float *,
               float *, int);
void   esvdcsint(double *, double *, void *, int, int *, double *,
              double *, int);
 
#ifdef __ESVERR
int esvscsin2_er( float *,  float *, void *, int, int, int,  float *,
                  float *, int, int, void *, int,  float *, int *);
int esvdcsin2_er(double *, double *, void *, int, int, int, double *,
                 double *, int, int, void *, int, double *, int *);
#else
void esvscsin2( float *,  float *, void *, int, int, int,  float *,
                float *, int, int, void *, int,  float *, int);
void esvdcsin2(double *, double *, void *, int, int, int, double *,
               double *, int, int, void *, int, double *, int);
#endif
 
 
/*  Numerical Quadrature Subroutines  */
 
void   esvsptnq( float *,  float *, int,  float *,  float *);
void   esvdptnq(double *, double *, int, double *, double *);
 
float  esvsglnq(void (*)(float*, float*, int*),  float,  float, int);
double esvdglnq(void (*)(float*, float*, int*), double, double, int);
 
float  esvsglnq2(void (*)(float*, int*, float*, int*, float*, int*), float,
                 float, int, float, float, int, void *, int);
double esvdglnq2(void (*)(float*, int*, float*, int*, float*, int*), double,
                 double, int, double, double, int, void *, int);
 
float  esvsglgq(void (*)(float*, float*, int*),  float,  float, int);
double esvdglgq(void (*)(float*, float*, int*), double, double, int);
 
float  esvsgraq(void (*)(float*, float*, int*),  float,  float, int);
double esvdgraq(void (*)(float*, float*, int*), double, double, int);
 
float  esvsghmq(void (*)(float*, float*, int*),  float,  float, int);
double esvdghmq(void (*)(float*, float*, int*), double, double, int);
 
 
/*  Random Number Generation Subroutines  */
 
void   esvsurand(double *, int,  float *);
void   esvdurand(double *, int, double *);
 
#ifdef __ESVERR
int esvsnrand_er(double *, int,  float *,  float *, int *);
int esvdnrand_er(double *, int, double *, double *, int *);
#else
void esvsnrand(double *, int,  float *,  float *, int);
void esvdnrand(double *, int, double *, double *, int);
#endif
 
void esvsurxor(int *, int,  float *,  float *);
void esvdurxor(int *, int, double *, double *);
 
/*  Parallel Processing Subroutines  */
 
void   esvdgemlp(void *, int, char *, void *, int, char *, void *, int, int,
                 int, int);
 
#ifdef __ESVERR
int esvdgefp_er(void *, int, int, int *, double *, int *);
 
int esvdppfp_er(double *, int, double *, int *);
 
int esvdskfsp_er(int, double *, int, int *, int *, double *, double *,
                 int *, void *, int, int);
 
int esvdgkfsp_er(int, double *, int, int *, double *, int, int *, int *,
                 double *, double *, int *, void *, int, int);
 
int esvscftp_er(int, void *, int, int, void *, int, int, int *,
                int, int, float, double *, int *, double *, int *);
 
int esvscft2p_er(int, void *, int, int, void *, int, int, int *, int *,
                 int, float, double *, int *, double *, int *);
 
int esvscft3p_er(void *, int, int, void *, int, int, int *, int *, int *,
                 int, float, double *, int *);
#else
void esvdgefp(void *, int, int, int *, double *, int);
void esvdppfp(double *, int, double *, int);
 
void esvdskfsp(int, double *, int, int *, int *, double *, double *,
               int, void *, int, int);
 
void esvdgkfsp(int, double *, int, int *, double *, int, int *, int *,
               double *, double *, int, void *, int, int);
 
void esvscftp(int, void *, int, int, void *, int, int, int,
              int, int, float, double *, int, double *, int);
 
void esvscft2p(int, void *, int, int, void *, int, int, int, int,
               int, float, double *, int, double *, int);
 
void esvscft3p(void *, int, int, void *, int, int, int, int, int,
               int, float, double *, int);
#endif
 
/*  Utility Subroutines  */
 
void   esvdsrsm(int, double *, int *, int *, int, int *, void *, void *, int);
void   esvstride(int, int, int *, char *, int);
void   esvivsset(int);
void   esveinfo(int, int *, int *);
int iessl(void);
 
#ifdef __ESVERR
int esvdgktrn_er(int, double *, int, int *, double *, int, int *, int *,
                 double *, int *);
int esvdsktrn_er(int, double *, int, int *, int *, double *, int *);
#else
void esvdgktrn(int, double *, int, int *, double *, int, int *, int *,
               double *, int);
void esvdsktrn(int, double *, int, int *, int *, double *, int);
#endif
void esvievops(int);
 
 

 

 
 

 
#endif
