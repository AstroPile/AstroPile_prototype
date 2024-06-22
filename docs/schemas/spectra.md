# Spectral information schema

| Field 	| Unit 	| Description 	|
|-------	|------	|-------------	|
| flux      | erg/cm^2/Angstrom/s   	| Flux            	|
| lambda    | Angstrom | Wavelengths   |
| ivar    	|      	| Inverse variance of the flux |
| lsf_sigma |       | Line spread function sigma / resolution |

## Normalization

By default, we should serve the flux as **not continuum-normalized**. We should add a post-processing function for normalization or have the normalized flux in a separate, extra field.

## Radial velocity correction

**Galaxy spextra** are often uncorrected. **Stellar spectra** are often corrected - let's serve them in this corrected form, but retain the **radial velocity** as an extra field.

## Additional parameters

Many downstream tasks may make use of additional, inferred parameters that are available in various surveys, such as inferred effective temperature, log g, mass, etc.