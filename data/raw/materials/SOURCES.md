# Raw material source manifest

## Au_Johnson.yml

- Retrieval URL: https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/master/database/data/main/Au/nk/Johnson.yml
- Retrieved on Windows execution host: 2026-09-13 (Asia/Shanghai)
- Upstream license statement: public domain, CC0 1.0 waiver in the downloaded file
- Original reference: P. B. Johnson and R. W. Christy, *Optical Constants of the Noble Metals*, Physical Review B 6, 4370-4379 (1972), DOI 10.1103/PhysRevB.6.4370
- State recorded upstream: room temperature
- Columns: vacuum wavelength in micrometres, refractive-index real part `n`, extinction coefficient `k`
- Interpolation rule for P0-P2: independent linear interpolation of `n` and `k`; no extrapolation
- SHA-256: populated and independently verified in run logs

## N-BK7

- Manufacturer source: SCHOTT N-BK7 optical glass datasheet
- Formula: SCHOTT Sellmeier equation with wavelength in micrometres
- Coefficients in `scripts/calspr/optics.py`: B1=1.03961212, B2=0.231792344, B3=1.01046945, C1=0.00600069867, C2=0.0200179144, C3=103.560653
- P0-P2 treatment: real-valued catalog dispersion; absorption and lot/temperature uncertainty are not modeled in the first pilot

## Connectivity record

The Kylin host attempted the Au retrieval URL directly on 2026-09-13 with `curl -L --max-time 20`. It timed out with curl error 28 after 20.001 s and zero received bytes. The same verified URL was therefore acquired on the Windows host and transferred by SCP. No mirror or policy bypass was used.

## Rakic-LD.yml

- Retrieval URL: https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/master/database/data/main/Au/nk/Rakic-LD.yml
- Original reference: A. D. Rakic et al., *Optical properties of metallic films for vertical-cavity optoelectronic devices*, Applied Optics 37, 5271-5283 (1998), DOI 10.1364/AO.37.005271
- Upstream description: Lorentz-Drude fit of experimental data from several sources
- Windows-retrieved SHA-256: `a91fca57ee99efc77354908fa11b9c4822e8852394042069caaac7610c52e602`
- Role: model-form/material-source sensitivity analysis; not an independent experiment performed in this project

## Olmon-ev.yml

- Retrieval URL: https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/master/database/data/main/Au/nk/Olmon-ev.yml
- Original reference: R. L. Olmon et al., *Optical dielectric function of gold*, Physical Review B 86, 235147 (2012), DOI 10.1103/PhysRevB.86.235147
- Upstream sample state: evaporated gold
- Windows-retrieved SHA-256: `be778621e6491fc4e2db6eee400fb329d44ada2b987f7bfd2e219c83fe32a338`
- Role: measured-film material sensitivity analysis; not a stand-in for the unknown device actually fabricated in any future experiment
- Acquisition note: the first local request was reset after about 21 s and produced no file; one bounded retry to the same official URL succeeded. The failed request is retained in the command history.
