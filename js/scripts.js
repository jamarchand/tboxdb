$(document).ready(function() {
            $('#example').DataTable();
            $('#example_filter').hide(); // Hide default search datatables where example is the ID of table
                     
                     
            //Set default values and initalize form
            $('input[id^=c]').hide();
            var index = $('#search_by_0').val().replace('B','');
            $('#colsearch_'+index+'_0').show();
            
                     
            //Add event listener to new row
            var field_name = '#search_by_0'
            $(field_name).on("change", function(){
                                      var rcount = document.getElementById("search_table").rows.length;
                                      $('input[id^=c]').hide();
                                      //Show only a select amount of fields
                                      
                                      for (i = 0; i < rcount; i++) {
                                      var index = $('#search_by_'+i.toString()).val().replace('B','');
                                      $('#colsearch_'+index+'_'+i.toString()).show();
                                      }
                                      
              });
            //$('#example').wrap('<div id="hide" style="display:none"/>');

        });
    
    
    
    function adv_search(){
        //Search
        var rcount = document.getElementById("search_table").rows.length;
        for (i = 0; i < rcount; i++) {
            var index = $('#search_by_'+i.toString()).val().replace('B','');
            //document.write(index)
            var col_val = $('#colsearch_'+index+'_'+i.toString()).val();

            $('#example')
                  .DataTable()
                  .columns(index)
                  .search(col_val);
        }
        
        $('#example')
            .DataTable()
            .draw();
            
            $('#search_visibility').hide();
            $('#table_visibility').show();
            
    }


           


//Set default values and initalize form
$('input[id^=c]').hide();
var index = $('#search_by_0').val().replace('B','');
$('#colsearch_'+index+'_0').show();

    
 //Adding new search rows
  function add_fields() {
      //Count number of rows
      var rcount = document.getElementById("search_table").rows.length;

      //Add fields
      document.getElementById("search_table").insertRow(rcount).innerHTML =  '<tr> <td class="pt-2"> <select class="custom-select" id="search_by_'+(rcount).toString()+'"> <option value="B0">TBID</option> <option value="B8">Sequence</option> <option value="B21">T-box class</option> <option value="B22">Regulation</option> <option value="B1">Host species</option> <option value="B2">Genomic accession</option> <option value="B3">Specifier Region</option> <option value="B4">Specifier</option> <option value="B5">T-box sequence</option> <option value="B20">Amino acid</option> <option value="B6">tRNA family</option> <option value="B7">Protein description</option> <option value="B13">Protein ID</option> <option value="B14">Host taxID</option> <option value="B15">Phylum</option> <option value="B16">Class</option> <option value="B17">Order</option> <option value="B18">Family</option> <option value="B19">Genus</option> </select> </td> <td class="pt-2"> <div class="pl-2"> <input id="colsearch_0_'+(rcount).toString()+'" placeholder="TBID (e.g. 0MHW7AY)" class="form-control"/> <input id="colsearch_1_'+(rcount).toString()+'" placeholder="Host organism (e.g. Bacillus cereus)" class="form-control"/> <input id="colsearch_2_'+(rcount).toString()+'" placeholder="Genomic accession of host organism (e.g. CP001176)" class="form-control" /> <input id="colsearch_3_'+(rcount).toString()+'" placeholder="Specifier region (e.g. AUCCA)" class="form-control" /> <input id="colsearch_4_'+(rcount).toString()+'" placeholder="Specifier (e.g. UCC)" class="form-control"/> <input id="colsearch_5_'+(rcount).toString()+'" placeholder="T-box sequence (e.g. UGGC)" class="form-control"/> <input id="colsearch_6_'+(rcount).toString()+'" placeholder="tRNA family (e.g. SER (GGA))" class="form-control"/> <input id="colsearch_7_'+(rcount).toString()+'" placeholder="Downstream protein description (e.g. serine--tRNA ligase)" class="form-control"/> <input id="colsearch_8_'+(rcount).toString()+'" placeholder="Sequence (e.g. CAGCTTATTCCCTGTCCTAGAAAGCCGGGGGTTGATGGAACCCGGTACACAATAAGCGAGCGAAATACACTCTGGAGTATCTCGGGTAATGCCGAGCGGAATCGCGAACGATAACTCGCTACGAGTGGCATGAATGGCCTGCGCAAGCGCTTGGGTCCATTGGTGCAAGATGGGTGGTAACACGGTTATTCAATCGTCCCTATGTCTACAATAGACAAGGGGCGATTTTTTGTTTTTTACCATAGGAGATTGAAGGAGTGTATGGCAGTGAATATTATCGACGAACTCGAATGGCGCGAAGCCGTCAATCA)" class="form-control"/> <input id="colsearch_13_'+(rcount).toString()+'" placeholder="Downstream protein ID (e.g. ACK61623)" class="form-control"/> <input id="colsearch_14_'+(rcount).toString()+'" placeholder="Host TaxID (e.g. 405532)" class="form-control"/> <input id="colsearch_15_'+(rcount).toString()+'" placeholder="Phylum (e.g. Firmicutes)" class="form-control"/> <input id="colsearch_16_'+(rcount).toString()+'" placeholder="Class (e.g. Bacilli)" class="form-control"/> <input id="colsearch_17_'+(rcount).toString()+'" placeholder="Order (e.g. Bacillales)" class="form-control"/> <input id="colsearch_18_'+(rcount).toString()+'" placeholder="Family (e.g. Bacillaceae)" class="form-control"/> <input id="colsearch_19_'+(rcount).toString()+'" placeholder="Genus (e.g. Bacillus)" class="form-control"/> <!--<input id="colsearch_20_'+(rcount).toString()+'" placeholder="tRNA sequence (e.g. GGTCCCGTGGTGTAGTGGTtAACATGCCTGCCTGTCACGCAGGAGAtCGCCGGTTCGACCCCGGTCGGGACCGCCA)" class="form-control"/>--> <input id="colsearch_20_'+(rcount).toString()+'" placeholder="Amino acid (e.g. Leu)" class="form-control"/> <input id="colsearch_21_'+(rcount).toString()+'" placeholder="T-box class (e.g. Class II)" class="form-control"/> <input id="colsearch_22_'+(rcount).toString()+'" placeholder="T-box regulation (e.g. Transcriptional)" class="form-control"/> </div> </td> </tr>';
        
    
    // document.getElementById("search_table").insertRow(1).innerHTML = '<tr><td class="pt-2"><select class="custom-select" id="search_by_'+(rcount).toString()+'"><option value="B0">TBDB_ID</option><option value="B1">Organism</option><option value="B2">Genomic Accession</option><option value="B4">Specifier Sequence</option></select></td><td class="pt-2"><div class="pl-2"><input id="colsearch_0_'+(rcount).toString()+'" placeholder="TBID" class="form-control"/><input id="colsearch_1_'+(rcount).toString()+'" placeholder="Organism" class="form-control"/><input id="colsearch_2_'+(rcount).toString()+'" placeholder="Genomic Accession" class="form-control"/><input id="colsearch_4_'+(rcount).toString()+'" placeholder="Specifier" class="form-control"/></div></td ></tr>';

     
      //Hide all fields
      $('input[id^=c]').hide();
      
      //Show only a select amount of fields
      for (i = 0; i <= rcount; i++) {
          var index = $('#search_by_'+i.toString()).val().replace('B','');
          $('#colsearch_'+index+'_'+i.toString()).show();
      }
      
      //Add event listener to new row
      var field_name = '#search_by_'+(rcount).toString()
      $(field_name).on("change", function(){
                                var rcount = document.getElementById("search_table").rows.length;
                                $('input[id^=c]').hide();
                                //Show only a select amount of fields
                                
                                for (i = 0; i <= rcount; i++) {
                                var index = $('#search_by_'+i.toString()).val().replace('B','');
                                $('#colsearch_'+index+'_'+i.toString()).show();
                                }
                                
        });
    }
