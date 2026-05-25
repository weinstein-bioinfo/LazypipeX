/* START JS ANNOTATION REPORT FILTERING */
	/* Datamodel for discrete and numetic filters */
	const dmodel = {
		tableid: NaN,
		table: null,
		coli: null,
		col_names: [],		// array of column names
		col_visib: null,		// column visibility values. Usage col_visib.get('colname')
		flt_disc: null,		// discrete filters as Map of Sets. Usage: flt_disc.get('dbtype').has(value)
		flt_num: null		// numeric filters as Map of numbers. Usage: flt_num.get('bitscore')
	};
		
	function annot_filters_init(table_id, prop_disc, prop_num){
		dmodel.tableid 			= table_id;
		dmodel.table 			= document.getElementById(table_id);
		dmodel.coli				= new Map();
		dmodel.col_names 		= [];
		dmodel.col_visib		 	= new Map();
		dmodel.flt_disc 			= new Map();
		dmodel.flt_num 			= new Map();
		
		const ths 		= dmodel.table.rows[0].getElementsByTagName("TH");
		
		// Reset coli, col_names, col_visib
		for( let coli=0; coli<ths.length; coli++){
			// datamodel
			let colname		= ths[coli].getAttribute("colname");
			dmodel.coli.set(colname,coli);
			dmodel.col_names.push(colname);
			dmodel.col_visib.set(colname, true);
			// elements
			let el = document.getElementById(colname+"_selector_a");
			if(el !== null){
				el.classList.add("selected");
			}
		}
		// Reset flt_disc
		for(let prop of prop_disc){
			let tmp = new Set();
			dmodel.flt_disc.set(prop,tmp);
		}
		const rows 		= dmodel.table.rows;
		for(i = 1; i<rows.length; i++){
			tds 		= rows[i].getElementsByTagName("TD");
			for(let prop of prop_disc){
				td		= tds[dmodel.coli.get(prop)];
				value 	= td.textContent || td.innerText;
				dmodel.flt_disc.get(prop).add(value);
			}
		}	// elements
		for(let prop of dmodel.flt_disc.keys()){
			for(let opt of dmodel.flt_disc.get(prop).values()){
				let el = document.getElementById(prop+'_'+opt+'_a');
				if(el !== null){
					el.classList.add("selected");
				}
			}
		}
		// Reset flt_num
		for(let prop of prop_num){
			dmodel.flt_num.set(prop,0.0);
			let el = document.getElementById(prop+'_slider');
			if(el !== null){
				el.value = 0.0;
			}
		}
		annot_filters_apply();
	}
	
	function annot_filters_reset(){
		// recall init with current table
		annot_filters_init(dmodel.tableid, Array.from(dmodel.flt_disc.keys()), Array.from(dmodel.flt_num.keys()) );
	}
	
	/* Sets value of a numeric filter */
	function set_filter_numeric(prop, value){
		var num = Number(value);
		if( dmodel.flt_num.has(prop) && !isNaN(num)){
			dmodel.flt_num.set(prop, num);	
			return true;
		}
		return false;
	}
	
	/* Sets values of a numeric filter and apply filtering */
	function preset_filter_numeric(value_map){
		if(value_map === null || value_map.constructor.name !== 'Map'){
			console.log("preset_num_values(): unexpected arg: "+value_map.constructor.name);
			return false;
		}
		for(let prop of value_map.keys()){
			if(dmodel.flt_num.has(prop)){
				let value 	= Number(value_map.get(prop));
				dmodel.flt_num.set(prop,value);
				let el 	= document.getElementById(prop+'_slider');
				if(el !== null){
					el.value = value;
				}
				let el2 = el.parentElement.nextElementSibling;
				if(el2 !== null){
					el2.value = value;
				}
			}
		}
		annot_filters_apply();
		return true;
	}
	
	/* Toggles value of a discrete filter to on/off */
	function toggle_filter_discrete(prop,value,el=null){
		if( dmodel.flt_disc.has(prop) ){
			dmodel.flt_disc.get(prop).has(value) ? dmodel.flt_disc.get(prop).delete(value) : dmodel.flt_disc.get(prop).add(value);
			if(el !== null){
				dmodel.flt_disc.get(prop).has(value) ? el.classList.add("selected") : el.classList.remove("selected");	
			}
			annot_filters_apply();
			return true;
		}
		return false;
	}

	
	function toggle_filter_column(prop,el=null){
		if(dmodel.col_visib.has(prop) ){
			dmodel.col_visib.get(prop) ? dmodel.col_visib.set(prop,false) : dmodel.col_visib.set(prop,true);
			if(el !== null){
				dmodel.col_visib.get(prop) ? el.classList.add("selected") : el.classList.remove("selected");
			}
			annot_filters_apply();
			return true;
		}
		return false;
	}
	
	function annot_filters_apply(){
		const rows = dmodel.table.rows;
		var i, tds, td, value;
		loop1: for(i = 1; i<rows.length; i++){
			tds 		= rows[i].getElementsByTagName("TD");
			for(let prop of dmodel.flt_disc.keys()){
				td		= tds[dmodel.coli.get(prop)];
				value 	= td.textContent || td.innerText;
				if( dmodel.flt_disc.get(prop).has(value) ){
					rows[i].style.display = "";
				}else{
					rows[i].style.display = "none";
					continue loop1;
				}
			}
			for(let prop of dmodel.flt_num.keys()){
				td		= tds[dmodel.coli.get(prop)];
				value 	= td.textContent || td.innerText;
				value	= Number(value);
				if( dmodel.flt_num.get(prop) <= value){
					rows[i].style.display = "";
				}else{
					rows[i].style.display = "none";
					continue loop1;
				}
			}
		}
		for(let j=0; j<dmodel.col_names.length; j++){
			let col = dmodel.col_names[j];
			if(dmodel.col_visib.get(col)){
				dmodel.table.querySelectorAll('tr th:nth-child('+(j+1)+')').forEach(el=>el.style.display = '');
				dmodel.table.querySelectorAll('tr td:nth-child('+(j+1)+')').forEach(el=>el.style.display = '');
			}
			else{
				dmodel.table.querySelectorAll('tr th:nth-child('+(j+1)+')').forEach(el=>el.style.display = 'none');
				dmodel.table.querySelectorAll('tr td:nth-child('+(j+1)+')').forEach(el=>el.style.display = 'none');
			}
		}
	}
/* END JS ANNOTATION REPORT NAVIGATION  */

/* START JS for NAVIGATION WITH SELECTABLE VIEWS */
/*
 * USAGE:
 *  HTML:
 *	<div class='reportview' id='view_id'> displayed at a time
 *	<p class='reportview_label'> navigation text updated with each displayed view
 *	<a class='reportview' id='view_id_a'> link highlighted on selection
 *
 *  JS:
 *  select_view(reportview_id,reportview_label)
 *  		- displays <div class='reportview'> with id $reportview_id, hides other reportviews.
 *  		- updates text in <p class='reportview_label'> to $reportview_label
*/ 
	let SELECTED_VIEW_ID		= "";
	let SELECTED_VIEW_LABEL	= "";
	
	function select_view(new_id,new_label){	
		let divs = document.getElementsByClassName("reportview");
		for(const div of divs){
			div.style.display="none";
		}
		if(document.getElementById( new_id ) !== null){
			document.getElementById( new_id ).style.display = "block";
		}
		if(document.getElementById(new_id+"_a")!== null){
			document.getElementById(new_id+"_a").className = "selected";
		}
		if(document.getElementById("selected_view_label") !== null){
			document.getElementById("selected_view_label").innerHTML = new_label;
		}
		
		if(document.getElementById( SELECTED_VIEW_ID + "_a" ) !== null){
			document.getElementById( SELECTED_VIEW_ID + "_a" ).className= "";
		}
		SELECTED_VIEW_ID 	= new_id;
		SELECTED_VIEW_LABEL	= new_label;
	}
/* END JS for NAVIGATION WITH SELECTABLE VIEWS */

/* START JS FOR TABLE SORTABLE */
function init_sortable_tables(){
	const tables = document.getElementsByTagName('TABLE');
	for (let i = 0; i < tables.length; i++) {
		// add onclick() listeners for sorting
		if(tables[i].className == "sortable"){
			const ths = tables[i].rows[0].getElementsByTagName("TH");
			if(!tables[i].id){ //make sure table has id
				tables[i].id = "table"+i;
			}
			for( let col=0; col<ths.length; col++){
				ths[col].onclick = function() { sortTable(tables[i].id,col) };
			}
		}
		// add classes to td-elements based on th types
		const ths 	= tables[i].rows[0].getElementsByTagName("TH");
		const types	= [];
		for( let col=0; col<ths.length; col++){
			types[col]	= ths[col].getAttribute('type');
		}
		const rows	= tables[i].rows;
		for(let i = 1; i<rows.length; i++){
			const tds = rows[i].getElementsByTagName("TD");
			for(let col=0; col<tds.length; col++){
				if(types[col]){
					tds[col].className = types[col];
				}
			}
		}
		/* sort table by cols marked as asc/desc
		for( let col=(ths.length-1); col>=0; col--){
			if(ths[col].className == "asc" || ths[col].className == "desc"){
				
			}
		}*/
	}
}
function sortTable(tableid,col) {
	var table 		= document.getElementById(tableid);
	var rows 		= table.rows;
	var dir 			= rows[0].getElementsByTagName("TH")[col].className;
	var col_type 	= rows[0].getElementsByTagName("TH")[col].getAttribute("type");
	if(dir == "asc"){        
		dir = "desc";
	}
	else if(dir=="desc"){
		dir = "asc";
	}
	else{
		if(col_type == "int" || col_type=="double" || col_type=="num"){
			dir = "desc";
		}
		else if(col_type =="string"){
			dir = "asc";
		}
		else{
			dir = "asc";
		}
	}
	/* Set sorting errows by setting class */
	var th_list = rows[0].getElementsByTagName("TH");
	for(let i =0; i< th_list.length; i++){
		th_list[i].className = "mixed";
	}
	th_list[col].className = dir;
	
	
    var i, x, y, shouldSwitch, switchcount = 0;
    var switching = true;
	while (switching) {
		// Start by saying: no switching is done:
		switching = false;
		for (i = 1; i < (rows.length - 1); i++) {
			shouldSwitch = false;
			
			if(col_type == "int" || col_type=="double" || col_type=="num"){
				x 	= rows[i].getElementsByTagName("TD")[col].innerHTML;
				y 	= rows[i + 1].getElementsByTagName("TD")[col].innerHTML;           
				x	= x.replace(/\s+/g,"");
				y	= y.replace(/\s+/g,"");
                tmp = parseFloat(x);
				if(!isNaN(tmp)){ x = tmp; }
				tmp = parseFloat(y);
				if(!isNaN(tmp)){ y = tmp; }
				if(x == "NA" || x == ""){ x= -Infinity; }
				if(y == "NA" || y == ""){ y= -Infinity; }				
			}
			else if(col_type == "element"){
				x 	= rows[i].getElementsByTagName("TD")[col].innerHTML;
				y 	= rows[i + 1].getElementsByTagName("TD")[col].innerHTML;
				const el_x = rows[i].getElementsByTagName("TD")[col].firstElementChild;
				const el_y = rows[i + 1].getElementsByTagName("TD")[col].firstElementChild;
				if( el_x){
					x	= el_x.innerHTML;
				}
				if( el_y){
					y	= el_y.innerHTML;
				}
			}
			else{
				x = rows[i].getElementsByTagName("TD")[col].innerHTML;
				y = rows[i + 1].getElementsByTagName("TD")[col].innerHTML;  			
			}
			if (dir == "asc") {
				if (x > y ) {
					shouldSwitch = true;
					break;
				}
			}
			else if (dir == "desc") {
				if (x < y) {
					shouldSwitch = true;
					break;}
			}
		}
		if (shouldSwitch) {
			rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
			switching = true;
			switchcount ++; 
		}    
	}
}
/* END JS FOR TABLE SORTABLE */
