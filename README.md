Code for:
[Synergy between CaMKII Substrates and beta-Adrenergic Signaling in Regulation of Cardiac Myocyte Ca(2+) Handling. Soltis AR, Saucerman JJ. Biophys J. 2010 Oct 6;99(7):2038-47. PMID: 20923637](https://www.cell.com/biophysj/fulltext/S0006-3495(10)00985-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0006349510009859%3Fshowall%3Dtrue)

<img src="images/SoltisBiophysJ2010_Figure1.jpg">

Contents:
- Main driver script: soltis_biophysJ2010_masterCompute.m
- ODE files:
    - soltis_biophysJ2010_masterODEfile.m - integrator script for 4 ODE model components
    - soltis_biophysJ2010_CaMKII_signaling_ODEfile.m - CaMKII signaling pathway ODEs
    - soltis_biophysJ2010_BARsignalling_odefile.m - Beta adrenergic signaling pathway ODEs
    - soltis_biophysJ2010_camODEfile.m - Ca/Calmodulin/CaMKII activation ODEs
    - soltis_biophysJ2010_eccODEfile.m - Excitation-contraction coupling (ECC) ODEs

Other:
    - initial_conditions - directory with initial condition states
    - clipdata2.m - utility for pulling data points from figures with Matlab


