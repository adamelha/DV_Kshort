# DV_Kshort

## The latest source code file is: Kshort_DDL_Niv_1Nov.cxx
## Kshort_DDL.cxx is the older file, which have many things but most of those are unecessary. Kshort_DDL_Niv_1Nov.cxx is the first iteration of Kshort_DDL.cxx
## Kshort_DDL.h is the header file to the source file. This file has been modified after its creation, but we don't need to worry about that one.
## jo_Kshort_DDL.py is the joboptions file created for Kshort analysis, although some portions have already been modified, still it needs to be modified
## Just to remind: the .cxx and .h files are placed in the src/ directory of your DDLStudies package whereas the .py should be present in share/ of DDLStudies package. 



####################### Steps suggested by Avner for Kshort analysis ######################################
1. Figure out how to add new variables to the ntuple created by the DDL code.
2. Write the Ks—> pi+ pi- reconstruction. A good way to do it is this:
        A. If the package already performs vertex fits to 2-tracks that aren’t leptons, use those vertices. If not, need to make the vertices by ourselves:
            a. Loop over all tracks not identified as leptons, index i 
            b. Loop over all other tracks, index j > i
            c. Perform a vertex fit on the 2 tracks. Figure out how to do that using standard ATLAS tools. # It seems that we won't need to do this anymore
        B. Require that the invariant mass (m_Ks) of the vertex is within, say, +/- 30 MeV of the Ks mass (from PDG). We will later tighten this cut, but for now this is good. The vertex is now a Ks candidate.
        C. Require that the Ks candidate vertex is displaced by at least 4 mm from the primary vertex (PV) in the radial direction 
        D. Calculate the vector {\vec d} from the PV to the Ks vertex and the angle alpha between \vec d and the Ks candidate momentum. Typically, a good cut to select a clean Ks sample is |alpha|<0.01. But after plotting this variable, you might select another cut.
3. Determine which radial regions to use. This can perhaps be understood from Jordi’s code. Otherwise, look at the code that implements the material veto. The regions to use are between material-veto regions.
4. Based on the m_Ks distribution, determine the final m_Ks cut (probably ~6 MeV, I’m guessing). After this last cut, the Ks sample should be very clean.
5. Count the number of Ks in each of the radial regions in the MC. Calculate r_MC(i) = (# of Ks in radial region i) / (# of Ks in radial region 0, which is 4 mm < r << the beampipe)
6. In the same way, calculate r_data(i)
7. Calculate R(i) = r_data(i) / r_MC(i)
8. How to use this R(i) with the killing factor - for that we need to look in the old documentation (from SUSY DV analysis).
   
#######################################################################################
