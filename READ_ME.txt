

--------------------------------------------------------------------------------------------------------

Description of the Supporting information, including scripts and data to generate the results and visualisations, from 
  'Temporary mate removal during incubation leads to variable compensation in a biparental shorebird'   
   by Bulla et al 2017
   https://doi.org/10.1101/117036

--------------------------------------------------------------------------------------------------------

WHEN USING any of the Supporting information, PLEASE CITE both the original paper and the Supporting information
	Bulla, M. at al. (2018).  Temporary mate removal during incubation leads to variable compensation in a biparental shorebird. bioRxive, https://doi.org/10.1101/117036. 
	Bulla, M. (2018).  Supporting information for 'Temporary mate removal during incubation leads to variable compensation in a biparental shorebird'. Open Science Framework, http://doi.org/10.17605/OSF.IO/mx82q ADD DATETIME.

For any publication making substantial use of the data, please contact Martin Bulla (bulla.mar@gmail.com), as the authors welcome the opportunity for collaboration and wish to comment prior to publication.

--------------------------------------------------------------------------------------------------------

CONTENT
1. Supplementary Figures, Tables, Methods and Ethics
2. Analyses
   2a. Scripts
   2b. Data
    a) csv files with metadata
    b) RData files prepared for analyses
3. Preprint

--------------------------------------------------------------------------------------------------------

1. Supplementary Figures, Tables, Methods and Ethics contain Supplementary items reffered to in the manuscript

--------------------------------------------------------------------------------------------------------

2a. Scripts (run with the below described datasets)
	- Functions&Constants.R necessary to run the  Analyses.R	
	- Analyses.R generates the resutls, figures and supporting information.
	
--------------------------------------------------------------------------------------------------------

2b. Data contain 
	
	a) coma delianated csv files # CHECK WHETHER ALL VARIABLES USED AT THE END
	   - captivity.csv - contains information about birds held in captivity
			author		: two letter abbreviation of the person taking care of the bird
			datetime	: datetime of the event
			box			: unique ID of the box, in which the bird was kept
			ring_num	: unique ring number of the bird
			phase		: phase of captivity condition: start: bird caught on the nest; middle - after ~12h; end - end of captivity period; released - bird was released in the vicinity of its nest
			fat			: fat scored on a qualitative scale (0-7) as the amount of fat in the axillary region and, if necessary, the furcular region (Meissner 2009)
			mass		: body mass of the bird in grams
			mass_worms	: mass in grams of meal worms given to the bird at start or left over at the end
			num_worms	: number of meal worms given to the bird at start or left over at the end
			mass_cat	: mass of cat food given to the bird at start or left over at the end
			sperm		: sperm taken using following sperm taking method: N - no sperm; Ms - abdominal massage and sperm taken from the feaces; S - sperm taken from feaces; M - abdominal massage, Clf - cloaca leavage  
			ringing_end	: datetime when ringing/measuring of the bird ended
			partner_present: was the bird's partner present during catching on the nest? Y - yes, N - no, other_bird - other birds present
			captivity_pk: captivity primary key
			remarks
		
	   - experiment_metadata.csv - contains summary of all dates and times related to the experiment
 		  nest				: unique identity of the nest
 		  IDfemale			: unique ID of the female
		  IDmale			: unique ID of the male
		  ID_taken			: ID of the bird that was taken into the captivity
		  ID_treated		: ID of the bird that was treated for compensation
 		  sex_taken			: sex of the bird that was taken into the captivity: f - female, m - male
 		  sex_treated		: sex of the bird that was treated for compensation: f - female, m - male
		  taken				: date and time when the bird was taken into captivity
		  ringing_end		: date and time when processing of the captured parent ended
		  start_			: date and time when the experimental period (control and treated) started (i.e. date time when experimental bird returned to the nest upon capture of its partner)
		  end_c				: date time indicating end of control experimental period 
		  end_				: date and time when the uniparental incubation of the treated parent ended (e.g. because its partner returned or because the treated parent deserted its nest) experimental period (control and treated) ended
		  release			: date and time when the captured parent was released from the captivity - demarkates end of experimental period
		  r_time			: indicates length of treated period in hours (that is difference between 'end_c' and 'release')
		  ### perhaps delete ### e_time			: indicates length of the treated period in hours (including also uniparental phase after release of the captive parent) - that is difference between 'end_' and 'end_c')
		  back				: has the released partner returned: y - yes, n-no
		  after				: what happened after release of the captive parent: the captive parent returned - return, or did not and its incubating partner continued incubating - continues or deserted - desertion
		  bout_a			: indicates how long it took the released captive parent to return to the nest in hours - i.e. difference between end_ and release (note that for plotting purposes for cases when the captive parent has not returned we hav assigned 16.84722222h)   
		  c_bout			: lenght of control incubation bout in minutes (etimated as median incubation bout of three pre-experimental bouts)
		  ####n
		  an				: indicates what happened if a parent has not returned
		  starved			: was the captive bird starved for the first 12h of captivity: yes, no
		  use_t				: indicates whether the nest was used for the analyses of compensation: y - yes, n - no
		  use_b				: indicates whether the nest was used for the analyses of the after experimental effects: y - yes, n - no
		  days				: indicates number of days a parent continued incubating (in case its partner has not returned
		  present			: was the bird's partner present during catching on the nest? Y - yes, N - no, other_bird - other birds present
		  lat				: latitude of the nest in decimals
		  lon				: longitude of the nest in decimals
		  end_state			: known end state of the nest: d - deserted; p - depredated; hg - hatching, hd - hatched, un - unknown, ud - undetermined, 1?hd,3d - 1 egg likely hatched, other three deserted,?p - probably depredated, ?hg - probably hatching
		  comments
		  mass_loss			: mass loss in gram whilein captivity
 		
	    - escape.csv 		- escape distances of experimental birds
		datetime_	: date and time of the visit to the nest when the escape distance was estimated
		nest		: identity of the nest 
		sex			: sex of the incubating parent:f-female, m-male
		bird_ID		: unique identification of the incubating parent
		type		: type of escape distance estimate: estimated - a person estimated its distance to the nest when the parent left; gps - a person marked its position when the parent left
		distm		: estimated distance of a person to the nest (in meteres) when the incubating parent left the nest
		
		 - escape_2011-2012.csv 		- escape distances of experimental birds
		year		: year of data collection
		datetime_	: date and time of the visit to the nest when the escape distance was estimated
		nest		: identity of the nest within given year
		sex			: sex of the incubating parent:f-female, m-male
		bird_ID		: unique identification of the incubating parent
		type		: type of escape distance estimate: estimated - a person estimated its distance to the nest when the parent left; gps - a person marked its position when the parent left
		distm		: estimated distance of a person to the nest (in meteres) when the incubating parent left the nest
		
	   - prop_inc.csv 		- contains mean and median proportion of incubation for a focal parent from three prior to treatment days with available data	
		nest		: identity of the nest of the focal parent
		mean_		: mean proportion of incubation
		med_		: median proportion of incubation
		

	   - incubation_start.csv - start of incubation for each nest
		year		: year of data collection
		nest		: identity of the nest within given year
		inc_start	: datetime when the incubation has started
		method		: method used to calculate the inc_start (see paper for details): laying - based on laying order, flotaion - estimated by floting eggs in the water

	   - bout_length_constancy_2011.csv
			nest		: identity of the nest 
			sex			: f-female, m-male of the incubating parent
			bout_start	: date and time when the particular incubation bout started
			bout_length	: length of particular incubation bout
			lat				: latitude of the nest in decimals
		    lon				: longitude of the nest in decimals
			
	   - non_experimental_nests.csv 	- fate of the nests that were not part of the experiment
		nest		: identity of the nest 
		end_state	: final known fate of the nest (? - indicates probable, numbers number of nests with given egg state): w - warm eggs, d - nest deserted, p - depredated, hs - start of hatching, hg - in process of hatching, hd - all chicks hatched, ln - chicks left the nest, un - unknown, ud - undetermined, b - egg has breaks (hatching will start soon), ho - eggs have holes (hatching is starting); bpd - a parent probably deserted, m - damaged or dried out egg
		comments	: additional information about 'state'

	   - cage.csv - indicates when cage was protecting a nest
		nest		: identity of the nest 
		on			: datetime when cage was put on the nest
		off			: datetime when cage was takenb off the nest

    c) Prepared for analyses 
	    - bout.Rdata - contains object 'b' with all incubation bouts of the studied nests 
		- experimental.Rdata - contains object 'b' (a subset of bout.Rdata) with only control and treatment bouts	
			"pk"             primary key - unique identity of each row
			"species"        sesa - indicates semipalmated sandpiper
			"nest"           unique identity of the nest
			"exper"         is the experimental bout control - c or treated t
			"bird_ID_filled" unique identity of the bird
			"sex"            genetic sex of the bird (f = female, m = male)
			"bout_ID"       unique ID of the incubation bout within nest
			"bout_type"      type of bout: inc - incubation bout, gap - exchange gap
			"inc_start"     date and time when the incubation at a given nest started
			"bout_start"    date and time when the given incubation bout or exchange gapstarted
			"bout_end"       date and time when the given incubation bout or exchange gap ended
			"bout_length"    length of incubation bout or exchange gap in minutes
			"inc_eff"     		incubation constancy based on temperature data   
			"inc_eff_2"     	incubation constancy based on RFID data   
			"inc_eff_3"      	incubation constancy based on temperature and RFID data   
			"disturb"        	how far in meters was the closest person (gps)
			"disturb_log"    ln(disturb)
			"t_ambient_med"		median tundra temperature in C (measure next to the nest)  during given bout
			"t_amb_avg"     mean tundra temperature in C based on loggers from the whole study are  during given bout
			"t_station_med" median ambient temperature in C (measure by the nearby weather stating)  during given bout
			"wind_sp_med"  median wind speed in m/s (measure by the nearby weather stating)  during given bout
			"wind_dir_med" direction of the wind in degrees (measure by the nearby weather stating)  during given bout
			"precip_sum"    precipitation in mm during given bout
			"h_med" 		median % humidity (measure next to the nest)  during given bout
			"h_avg_med"		mean % humidity based on loggers from the whole study are  during given bout
			"light_int_med"	median light intensity if availalbe
			"sun_elev_med"  sun elevation in degrees
			"t_nest_med"	median nest temperature in C  during given bout
			"t_type"		type of T probe in the nest
			"t_amb_type"    type of T probe next to the nest

	    - constancy_for_Fig_3.Rdata - contains object 'u' with data to generate Fig 3 
				pk    : primary key - unique row ID
				datetime_: date and time of the event
				inc_t        : temperature based incubation: 1 = incubation, 0 = no incubation
				incubation        : RFID based incubation: 1 = incubation, 0 = no incubation
				datetime_c			: indicates hours to (minus) or from the end of control experimental period
				roll_con        : 12min rolling constancy of incubation (proportion)
				nest		: identity of the nest 
				col_line and col_	: indicates color of the line in Figure 4
				n					: indicates order of the nest
				sex		: sex of the incubating parent:f-female, m-male

--------------------------------------------------------------------------------------------------------


	

