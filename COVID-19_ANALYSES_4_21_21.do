
********************************************************************************
//REPLICATION CODE FOR NON-CAUSAL COVID-19 NETWORK ANALYSES
//LAST UPDATED 4-21-21
********************************************************************************
clear 
clear matrix

*********************************************************************************
*********************************************************************************
//MAIN PANEL ANALYSES
*********************************************************************************
*********************************************************************************

//TABLE 1: DESCRIPTIVE STATS
//note that this needs to be manually coppied into excel - does not generate a matrix for export

//copied onto RESULTS.xls, sheet(descriptives)

use "${root}\covid_2_24_21_panel_cleaned_no_AK.dta", clear 

summarize cases deaths case_rate death_rate n_covid n_delta w_covid w_delta tl_covid disadvantage p65_over abavg_nhw abavg_nhb abavg_hisp ppopurban

********************************************************************************

//TABLES 2 AND 3: FIXED EFFECTS NBREG MODELS
//includes AIC, BIC

use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

cd "${root}"

xtset geoid time

foreach outcome in "cases" "deaths" {
	
	putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(`outcome'_FE) modify

	//space only
	xtnbreg `outcome' std_w_covid ib1.time, exposure(total_pop) fe
	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel A1 = matrix(result), names nformat(number_d2) hcenter

	estat ic
	matrix IC = r(S)
	putexcel A40 = matrix(IC)

	//network only
	xtnbreg `outcome' std_n_covid ib1.time, exposure(total_pop) fe
	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel I1 = matrix(result), names nformat(number_d2) hcenter

	estat ic
	matrix IC = r(S)
	putexcel I40 = matrix(IC)

	//space and network
	xtnbreg `outcome' std_n_covid std_w_covid ib1.time, exposure(total_pop) fe
	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel Q1 = matrix(result), names nformat(number_d2) hcenter

	estat ic
	matrix IC = r(S)
	putexcel Q40 = matrix(IC)

	//space, network, and all controls
	xtnbreg `outcome' std_n_covid std_w_covid std_n_delta std_w_delta std_tl_covid ib1.time, exposure(total_pop) fe
	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel Y1 = matrix(result), names nformat(number_d2) hcenter

	estat ic
	matrix IC = r(S)
	putexcel Y40 = matrix(IC)

}

********************************************************************************
//MULTILEVEL RANDOM INTERCEPT (MENBREG) MODELS
//includes AIC, BIC, MAAPE 
//note that estimates are used for coefficient plot figures (same excel sheets)
//note that these esimates are used for fit stat table

use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

cd "${root}"

foreach outcome in "cases" "deaths" {
		
	putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(`outcome') modify

	//space and controls
	menbreg `outcome' std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel A1 = matrix(result), names nformat(number_d2) hcenter

	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
	scalar MAAPE_`outcome'_1 = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time

	estat ic
	matrix IC = r(S)
	putexcel A63 = matrix(IC)
	putexcel A65 = ("MAAPE")
	putexcel B65 = (MAAPE_`outcome'_1), nformat(number_d2)

	//network and controls
	menbreg `outcome' std_n_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:

	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel I1 = matrix(result), names nformat(number_d2) hcenter

	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
	scalar MAAPE_`outcome'_2 = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time

	estat ic
	matrix IC = r(S)
	putexcel I63 = matrix(IC)
	putexcel I65 = ("MAAPE")
	putexcel J65 = (MAAPE_`outcome'_2), nformat(number_d2)

	//space, network, and controls
	menbreg `outcome' std_n_covid std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:

	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel Q1 = matrix(result), names nformat(number_d2) hcenter

	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
	scalar MAAPE_`outcome'_3 = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time

	estat ic
	matrix IC = r(S)
	putexcel Q63 = matrix(IC)
	putexcel Q65 = ("MAAPE")
	putexcel R65 = (MAAPE_`outcome'_3), nformat(number_d2)
}

********************************************************************************
// PERMUTATION TEST TABLE BASED ON FULL RANDOM INTERCEPT MODEL (TABLE 4)
//note that these are computationally intensive, so normally run on the PRI server

use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

//define program
program define MAAPEprog, rclass //purpose of this program is to generate MAE for a model, given an outcome
    version 16.1
    args outcome
	menbreg `outcome' std_n_covid std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
    return scalar MAAPE = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time
end


cd "${root}"

foreach outcome in "cases" "deaths" {

	putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(permute_`outcome') modify
	local excel_row = 2

	putexcel A1 = ("variable")
	putexcel B1 = ("proportion")
	putexcel C1 = ("trials")

	foreach permutev in std_n_covid std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban {
		permute `permutev' r(MAAPE), reps(1000) seed(4815) left: MAAPEprog "`outcome'"
				
		putexcel A`excel_row' = ("`permutev'")
		mat result = r(p)
		putexcel B`excel_row' = matrix(result), nformat(number_d2) hcenter
		mat result = r(reps)
		putexcel C`excel_row' = matrix(result), nformat(number_d2) hcenter
		local excel_row = `excel_row' + 1
				
	}
}

********************************************************************************
********************************************************************************
//MENBREG SUPPLEMENTAL TABLES WITH LEVEL-1 COVARIATES

use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

cd "${root}"

foreach outcome in "cases" "deaths" {
	
	putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(`outcome'_delta) modify

	menbreg `outcome' std_n_covid std_w_covid std_n_delta std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:

	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel A1 = matrix(result), names nformat(number_d2) hcenter

	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
	scalar MAAPE_`outcome' = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time

	estat ic
	matrix IC = r(S)
	putexcel A43 = matrix(IC)
	putexcel A45 = ("MAAPE")
	putexcel B45 = (MAAPE_`outcome'), nformat(number_d2)

}
*******************************************************************************
use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

tempfile panel
save `panel'

//spatiall lagged outcomes

cd "${SHAPE_root}"

//unzipfile tl_2015_us_county.zip
//spshape2dta tl_2015_us_county
//you only need these lines if all the shapefile stuff is in a zipped folder

use tl_2015_us_county.dta, clear
destring GEOID, gen(geoid)
//don't need this because already made it in a previous run

spmatrix create contiguity W, normalize(row) 

merge 1:m geoid using `panel'
drop if _merge == 2 //otherwise stuff gets messed up.
drop if _merge == 1
drop _merge

cd "${root}"

xtset _ID time

putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(spxt_reg) modify

spxtregress std_death_rate std_n_covid ib1.time, fe dvarlag(W) errorlag(W) force
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter

spxtregress std_case_rate std_n_covid ib1.time, fe dvarlag(W) errorlag(W) force
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

*******************************************************************************
//UPDATED 3-13-21

use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

cd "${root}"

foreach x of numlist 1 8 16 23 {
    
    use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 
    keep if time == `x'
	
	xtile n_covid_`x' = n_covid, nq(10)
	xtile death_rate_`x' = death_rate, nq(10)
	xtile case_rate_`x' = case_rate, nq(10)
	keep geoid n_covid_`x' death_rate_`x' case_rate_`x'
	
	tempfile deciles_`x'
	save `deciles_`x''
}

use `deciles_1', clear
foreach x of numlist 8 16 23 {
	merge 1:1 geoid using `deciles_`x''
	drop _merge
}

tempfile all
save `all'

cd "${SHAPE_root}" //adding back in the string version of geoid makes the merge into arcgis better

use tl_2015_us_county.dta, clear
destring GEOID, gen(geoid)
keep GEOID geoid 

merge 1:1 geoid using `all'
drop _merge

export excel using "${root}\MAIN_MAP_UPDATED.xls", firstrow(variables) replace
//this will then be used in ARCGIS (available on PRI server) to make MAIN figures 1, 2, and 5
//finished figures are in the MAIN_FIGURES_MAPS folder

********************************************************************************
//analyses by disadvantage thirtiles

use "${root}\covid_2_24_21_panel_cleaned_no_AK", clear 

cd "${root}"

gen dis_th = .
centile std_disadvantage, centile (33.3333 66.6666)
replace dis_th = 1 if std_disadvantage <= r(c_1)
replace dis_th = 2 if (std_disadvantage > r(c_1)) & (std_disadvantage <= r(c_2))
replace dis_th = 3 if (std_disadvantage > r(c_2)) & (std_disadvantage != .)

putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(low_dis) modify
menbreg deaths std_n_covid std_n_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time if dis_th == 1, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter
menbreg cases std_n_covid std_n_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time if dis_th == 1, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(mid_dis) modify
menbreg deaths std_n_covid std_n_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time if dis_th == 2, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter
menbreg cases std_n_covid std_n_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time if dis_th == 2, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(high_dis) modify
menbreg deaths std_n_covid std_n_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time if dis_th == 3, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter
menbreg cases std_n_covid std_n_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time if dis_th == 3, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

*********************************************************************************
*********************************************************************************
//SPLIT NETWORK ANALYSES (CONTIGUOUS AND NON-CONTIGUOUS NETWORK)
*********************************************************************************
*********************************************************************************

use "${root}\covid_2_24_panel_split_ntwk", clear 

cd "${root}"

foreach outcome in "cases" "deaths" {
		
	putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(`outcome'_split_ntwk) modify

	menbreg `outcome' std_n1_covid std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:

	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel A1 = matrix(result), names nformat(number_d2) hcenter

	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
	scalar MAAPE_`outcome' = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time

	estat ic
	matrix IC = r(S)
	putexcel A42 = matrix(IC)
	putexcel A44 = ("MAAPE")
	putexcel B44 = (MAAPE_`outcome'), nformat(number_d2)

	menbreg `outcome' std_n0_covid std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:

	mat result = r(table)
	matrix result = result[1..6,1...]'
	putexcel I1 = matrix(result), names nformat(number_d2) hcenter

	predict phat //predicted counts for given outcome
	gen resid = phat - `outcome' //creates residuals
	gen absresid = abs(resid) //creates absolute value of residuals
	gen rel_resid = atan(absresid/`outcome')
	replace rel_resid = 1.5707963 if rel_resid ==.
	replace rel_resid = . if e(sample) == 0
	egen sumresid = total(rel_resid) 
	qui sum sumresid //summarizes sumresid
	scalar sumres = r(mean) //pulls the mean value for sumresid, which I can define as a scalar
	scalar MAAPE_`outcome' = sumres/e(N) //divide sumresid by total N for the model, and return this value as MAE
	drop phat resid absresid rel_resid sumresid //clean up all varnames for next time
	scalar drop sumres //clean up scalar name for next time

	estat ic
	matrix IC = r(S)
	putexcel I42 = matrix(IC)
	putexcel I44 = ("MAAPE")
	putexcel J44 = (MAAPE_`outcome'), nformat(number_d2)

}

*********************************************************************************
*********************************************************************************
//STRATIFIED NETWORK ANALYSES (THREE INCOME LEVELS)
*********************************************************************************
*********************************************************************************

//predictive models with network stratified by income for supplemental tables

use "${root}\covid_2_24_panel_income", clear

cd "${root}"

**low income
putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(low_inc) modify

menbreg deaths std_n_1_covid std_n_1_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter

menbreg cases std_n_1_covid std_n_1_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

**middle income
putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(mid_inc) modify

menbreg deaths std_n_2_covid std_n_2_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter

menbreg cases std_n_2_covid std_n_2_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

**high income
putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(high_inc) modify

menbreg deaths std_n_3_covid std_n_3_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter

menbreg cases std_n_3_covid std_n_3_delta std_w_covid std_w_delta std_tl_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban ib1.time, exposure(total_pop) vce(cluster STATEFP) || geoid:
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

*********************************************************************************
*********************************************************************************
//CUMULATIVE DATA ANALYSES
*********************************************************************************
*********************************************************************************

//LOO-MAAPE for fit stat table
// note that this takes a while to run, previously done on PRI server

use "${root}\covid_2_24_21_cross_section_cumulative_no_AK.dta", clear

cd "${root}"

putexcel set covid_results_UPDATED_2_24_no_AK.xls, sheet(LOO_w_only) modify
local excel_row = 1
foreach outcome in "cases" "deaths" {
	local sum_error = 0
	use covid_2_24_21_cross_section_cumulative_no_AK.dta, clear
	local N = _N
	forval n = 1(1)`N' {
		drop if _n == `n'
		menbreg `outcome' std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban, exposure(total_pop) vce(robust) || STATEFP:	
		use covid_2_24_21_cross_section_cumulative_no_AK.dta, clear
		predict phat
		summarize phat if _n == `n'
		scalar phat_n = r(mean)
		assert phat_n !=.
	
		summarize `outcome' if _n == `n'
		scalar obs_n = r(mean)
		assert obs_n !=.
	
		scalar error_n = abs(phat_n - obs_n)/obs_n 
	
		if error_n==. {
			scalar error_n = 8e+307 //approximate infinitiy if obs_n = 0
		}
		
		assert obs_n !=0 | error_n == 8e+307
	
		local sum_error = `sum_error' + atan(error_n) 
	
		scalar drop phat_n obs_n error_n
	
		disp `sum_error'
		disp `N'
	}
	
	scalar MAAPE = `sum_error' / `N'
	putexcel A`excel_row' = ("`outcome'")
	putexcel B`excel_row' = (MAAPE), nformat(number_d2) hcenter
	local excel_row = `excel_row' + 1
}
**********************************************************************************

putexcel set covid_results_UPDATED_2_24_no_AK.xls, sheet(LOO_n_only) modify
local excel_row = 1
foreach outcome in "cases" "deaths" {
	local sum_error = 0
	use covid_2_24_21_cross_section_cumulative_no_AK.dta, clear
	local N = _N
	forval n = 1(1)`N' {
		drop if _n == `n'
		menbreg `outcome' std_n_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban, exposure(total_pop) vce(robust) || STATEFP:	
		use covid_2_24_21_cross_section_cumulative_no_AK.dta, clear
		predict phat
		summarize phat if _n == `n'
		scalar phat_n = r(mean)
		assert phat_n !=.
	
		summarize `outcome' if _n == `n'
		scalar obs_n = r(mean)
		assert obs_n !=.
	
		scalar error_n = abs(phat_n - obs_n)/obs_n 
	
		if error_n==. {
			scalar error_n = 8e+307 //approximate infinitiy if obs_n = 0
		}
		
		assert obs_n !=0 | error_n == 8e+307
	
		local sum_error = `sum_error' + atan(error_n) 
	
		scalar drop phat_n obs_n error_n
	
		disp `sum_error'
		disp `N'
	}
	
	scalar MAAPE = `sum_error' / `N'
	putexcel A`excel_row' = ("`outcome'")
	putexcel B`excel_row' = (MAAPE), nformat(number_d2) hcenter
	local excel_row = `excel_row' + 1
}
	
********************************************************************************

putexcel set covid_results_UPDATED_2_24_no_AK.xls, sheet(LOO_w_and_n) modify
local excel_row = 1
foreach outcome in "cases" "deaths" {
	local sum_error = 0
	use covid_2_24_21_cross_section_cumulative_no_AK.dta, clear
	local N = _N
	forval n = 1(1)`N' {
		drop if _n == `n'
		menbreg `outcome' std_n_covid std_w_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban, exposure(total_pop) vce(robust) || STATEFP:	
		use covid_2_24_21_cross_section_cumulative_no_AK.dta, clear
		predict phat
		summarize phat if _n == `n'
		scalar phat_n = r(mean)
		assert phat_n !=.
	
		summarize `outcome' if _n == `n'
		scalar obs_n = r(mean)
		assert obs_n !=.
	
		scalar error_n = abs(phat_n - obs_n)/obs_n 
	
		if error_n==. {
			scalar error_n = 8e+307 //approximate infinitiy if obs_n = 0
		}
		
		assert obs_n !=0 | error_n == 8e+307
	
		local sum_error = `sum_error' + atan(error_n) 
	
		scalar drop phat_n obs_n error_n
	
		disp `sum_error'
		disp `N'
	}
	
	scalar MAAPE = `sum_error' / `N'
	putexcel A`excel_row' = ("`outcome'")
	putexcel B`excel_row' = (MAAPE), nformat(number_d2) hcenter
	local excel_row = `excel_row' + 1
}
********************************************************************************
use "${root}\covid_2_24_21_cross_section_cumulative_no_AK.dta", clear

//spatial autoregressive models using cumulative data.
//for supplemental tables.
//note that spatial matrix W needs to be loaded ahead of time.

cd "${root}"

putexcel set RESULTS_UPDATED_2_24_no_AK.xls, sheet(sp_reg) modify

spregress std_case_rate std_n_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban i.STATEFP, ml dvarlag(W) errorlag(W) vce(robust) force
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel A1 = matrix(result), names nformat(number_d2) hcenter
spregress std_death_rate std_n_covid std_disadvantage std_p65_over abavg_nhw abavg_nhb abavg_hisp std_ppopurban i.STATEFP, ml dvarlag(W) errorlag(W) vce(robust) force
mat result = r(table)
matrix result = result[1..6,1...]'
putexcel I1 = matrix(result), names nformat(number_d2) hcenter

********************************************************************************

