#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <curl/curl.h>
#include <unistd.h>


// Function declarations
void print_help();
void output_annotation_header(FILE *oah_annotated_output_fp);  // Outputs column header for annotation file
void free_2d_str_array(char **f2sa_array, int f2sa_array_len);  // Frees memory for a 2-dimensional char array
int main(int argc, char *argv[]); 
// END OF: Function declarations



// global (begin with g_) strings and variables
char *g_version_name = "Tempus annotator, Version 1.0\n";
char g_dp_format_str[3] = "DP";  // VCF format field tag/id for "Read depth"
char g_dpr_format_str[4] = "DPR";  // VCF format field tag/id for "Number of observations for each allele"
char g_ro_format_str[3] = "RO";  // VCF format field tag/id for "Reference allele observation count"
char *g_url_address = "http://exac.hms.harvard.edu";  // URL address for EXAC database
char *g_rest_variant_variant = "/rest/variant/variant/";  // REST call for EXAC database; retrieves variant information
char *g_rest_variant_ordered_csqs = "/rest/variant/ordered_csqs/";  // REST call for EXAC database; retrieves ordered list of variant consequences (most severe to least severe)

int g_so_term_str_max = 1024;  // Maximum length for an SO string

char *g_help_str_array[10] = 
{
  "\nUsage: tempus -i <VCF-formatted input file> -o <Annotated output file> [optional parameters]\n",
  "\nRequired Parameters:\n",
  "\t-i \tVCF-formatted input file\n",
  "\t-o \tAnnotated output file\n",
  "\nOptional Parameters:\n",
  "\t-s \tSO terms file [SO_terms_sorted_by_ensembl_estimated_severity.txt]\n",
  "\nOther Parameters:\n",
  "\t-h \tPrint help\n",
  "\t-v \tPrint version\n",
  "\nSee README file for additional help.\n"
};
int g_help_str_array_len = sizeof(g_help_str_array)/sizeof(g_help_str_array[0]);

char *g_annotated_header_str_array[10] =
{
  "Chromosome",
  "\tPosition",
  "\tRef",
  "\tAlt",
  "\tDepth (VCF DP format field)",
  "\tReads supporting Alt (VCF DPR format field)",
  "\tPercent Reads supporting Alt vs. Reads supporting Ref (VCF)",
  "\tAllele Frequency (EXAC)",
  "\tVariant Type (VCF TYPE info field)",
  "\tVariant Type/Consequence (EXAC)"
};
int g_annotated_header_str_array_len = sizeof(g_annotated_header_str_array)/sizeof(g_annotated_header_str_array[0]);
// END OF: global (begin with g_) strings and variables


struct memory_struct 
{
  char *memory;
  size_t size;
};
 
static size_t


// Allocates and stores curl downloaded data    
curl_write_function(char *in_chunk, size_t in_size, size_t in_num, void *userp)
{
  size_t in_memory_size = in_size * in_num;
  struct memory_struct *mem = (struct memory_struct *)userp;
 
  mem->memory = realloc(mem->memory, mem->size + in_memory_size + 1);
  if(mem->memory == NULL)
  {
    printf("\nUnable to re-allocate memory for curl_chunk.\n");
    return 0;
  }
 
  memcpy(&(mem->memory[mem->size]), in_chunk, in_memory_size);
  mem->size += in_memory_size;
  mem->memory[mem->size] = 0;
 
  return in_memory_size;
}
 

void print_help()
{
  int ph_loop;
  
  printf("\n%s", g_version_name);
  for(ph_loop=0;ph_loop<g_help_str_array_len;ph_loop++)
  {
    printf("%s", g_help_str_array[ph_loop]);
  }
  printf("\n");
  return;
}


// Outputs column header for annotation file
void output_annotation_header(FILE *annotated_output_fp)
{
  int oah_loop;
  for(oah_loop=0;oah_loop<g_annotated_header_str_array_len;oah_loop++)
  {
    fprintf(annotated_output_fp, "%s", g_annotated_header_str_array[oah_loop]);
  }
  fprintf(annotated_output_fp, "\n");
  return;
}


// Frees memory for a 2-dimensional char array
void free_2d_str_array(char **f2sa_array, int f2sa_array_len)
{
  int f2sa_loop = 0;
  for(f2sa_loop=0;f2sa_loop<f2sa_array_len;f2sa_loop++)
  {
    free(f2sa_array[f2sa_loop]);
  }
  free(f2sa_array);
  return;
}


int main(int argc, char *argv[]) 
{
  
  // File handles and file names
  FILE *vcf_input_fp, *annotated_output_fp, *so_terms_fp;
  
  char *vcf_input_file_name = NULL;
  char *annotated_output_file_name = NULL;
  char *so_terms_file_name = NULL;
  // END OF: File handles and file names

  
  int m_a_loop = 0;  // Used in for loops
  int m_so_terms_count = 0;
  
  char vcf_chr_name[1024], vcf_id[1024], vcf_ref[1024], vcf_alt[1024], vcf_qual[1024], vcf_filter[1024], vcf_info[1024], vcf_format[1024], vcf_normal_sample[1024], vcf_sample[10000]; // assume 1 tumor sample, extra samples stored in rvif_sample but ignored
  int m_pos_start;  // Variant start position in reference
  
  char **m_so_terms, ** m_so_terms_no_quotes;
  
  int m_read_variants = 0;
  char m_temp_line[10000];
    
  char *m_pchr;
  int m_so_terms_len = 0;

  char url[1024];
  CURL *curl_handle;
  CURLcode res;
 
  struct memory_struct curl_chunk;
  

  // Read command line parameters
  printf("\n");  // line spacer after command line

  int opt = 0;
  while ((opt = getopt (argc, argv, "i:o:s:vh")) != -1) 
  {
    switch (opt)
    {
      case 'i':
	vcf_input_file_name = optarg;
	break;
      case 'o':
	annotated_output_file_name = optarg;
	break;
      case 's':
	so_terms_file_name = optarg;
	break;
      // v1_1
      case 'h':
	print_help();
	return 0;
      case 'v':
	printf("%s\n", g_version_name);
	return 0;

      case '?':
	print_help();
	if (  optopt == 'i' || optopt == 'o' || optopt == 's' )
	  fprintf (stderr, "\nERROR: Option -%c requires an argument.\n\n", optopt);
	else if (isprint (optopt))
	  fprintf (stderr, "\nERROR: Unknown option `-%c'.\n\n", optopt);
	else
	  fprintf (stderr, "\nERROR: Unknown option character `\\x%x'.\n\n", optopt);
	return 1;
      default:
	print_help();
	abort ();
    }
  }
  printf("\n%s\n\n", g_version_name);
  printf("VCF input file: %s\n", vcf_input_file_name);
  printf("Annotated output file: %s\n", annotated_output_file_name);
  if( so_terms_file_name == NULL )
  {
    so_terms_file_name = "SO_terms_sorted_by_ensembl_estimated_severity.txt";
  }
  printf("SO terms file: %s\n", so_terms_file_name);
  // END OF: Read command line parameters
        

  // prepare for curl request
  curl_chunk.memory = malloc(1);  // will be grown as needed in curl_write_function
  curl_chunk.size = 0;  // downloaded memory size
 
  curl_global_init(CURL_GLOBAL_ALL);
  curl_handle = curl_easy_init();  // initiate curl session
  curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, curl_write_function);  // set for data to be sent to curl_write_function
  curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *) &curl_chunk);  // pass curl_chunk structure to curl_write_function
  curl_easy_setopt(curl_handle, CURLOPT_USERAGENT, "libcurl-agent/1.0");  // user-agent field
  // END OF: prepare for curl request
 
 
  // Open and count lines in SO terms file
  if( so_terms_file_name != NULL )
  {
      so_terms_fp = fopen(so_terms_file_name, "r");
      if( so_terms_fp == NULL )
      {
	  printf("\nERROR: Could not open %s\n\n", so_terms_file_name);
	  return 1;
      }
      fgets(m_temp_line, sizeof(m_temp_line), so_terms_fp);  // ignore header
      while( fgets(m_temp_line, sizeof(m_temp_line), so_terms_fp) )
      {
	m_so_terms_len += 1;  // count lines
      }
  }
  else
  {
      print_help();
      printf("\nERROR: No SO terms file specified.\n\n");
      return 1;
  }
  // END OF: Open and count lines in SO terms file

  // Allocate m_so_terms
  m_so_terms = malloc(m_so_terms_len*sizeof(char*));
  if( m_so_terms == NULL ) 
  {
    printf("\nUnable to allocate memory for SO terms.\n");
    return 1;
  }
  for(m_a_loop=0;m_a_loop<m_so_terms_len;m_a_loop++)
  {
    m_so_terms[m_a_loop] = malloc((g_so_term_str_max+1)*sizeof(char));
    if( m_so_terms[m_a_loop] == NULL )  
    {
      printf("\nUnable to allocate memory for SO terms string\n");
      return 1;
    }
  }
  // END OF: Allocate m_so_terms

  // Allocate m_so_terms_no_quotes
  m_so_terms_no_quotes = malloc(m_so_terms_len*sizeof(char*));
  if( m_so_terms_no_quotes == NULL ) 
  {
    printf("\nUnable to allocate memory for SO terms (no quotes).\n");
    return 1;
  }
  for(m_a_loop=0;m_a_loop<m_so_terms_len;m_a_loop++)
  {
    m_so_terms_no_quotes[m_a_loop] = malloc((g_so_term_str_max+1)*sizeof(char));
    if( m_so_terms_no_quotes[m_a_loop] == NULL )  
    {
      printf("\nUnable to allocate memory for SO terms (no quottes) string\n");
      return 1;
    }
  }
  // END OF: Allocate m_so_terms_no_quotes

  // Read SO terms file
  m_so_terms_count = 0;
  rewind(so_terms_fp);
  fgets(m_temp_line, sizeof(m_temp_line), so_terms_fp);  // ignore header
  while( fgets(m_temp_line, sizeof(m_temp_line), so_terms_fp) )
  {
    m_pchr = strtok(m_temp_line, "\t");
    memcpy(m_so_terms_no_quotes[m_so_terms_count], m_pchr, strlen(m_pchr)+1);
    memcpy(m_so_terms[m_so_terms_count] + 1, m_so_terms_no_quotes[m_so_terms_count], strlen(m_pchr));
    m_so_terms[m_so_terms_count][0] = '"';
    m_so_terms[m_so_terms_count][strlen(m_pchr) + 1] = '"';
    m_so_terms[m_so_terms_count][strlen(m_pchr) + 2] = '\0';
    m_so_terms_count += 1;
  }
  printf("\nSO terms read: %d\n", m_so_terms_len);
  fclose(so_terms_fp);
  // END OF: Read SO terms file



  // Open annotated output file
  if( annotated_output_file_name != NULL )
  {
      annotated_output_fp = fopen(annotated_output_file_name, "w");

      if(annotated_output_fp == NULL)
      {
	  printf("\nERROR: Could not open %s\n\n", annotated_output_file_name);
	  return 1;
      }
  }
  else
  {
      print_help();
      printf("\nERROR: No output file specified.\n\n");
      return 1;
  }
  // END OF: Open annotated output file
      
		
  // Open and find length of variant list file
  if( vcf_input_file_name != NULL )
  {
      vcf_input_fp = fopen(vcf_input_file_name, "r");
      if( vcf_input_fp == NULL )
      {
	  printf("\nERROR: Could not open %s\n\n", vcf_input_file_name);
	  return 1;
      }
      // END OF: Estimate variants, tempus v1_0
  }
  else
  {
      print_help();
      printf("\nERROR: No VCF input file specified.\n\n");
      return 1;
  }
  // END OF: Open and find length of variant list file


  // Output annotation header
  output_annotation_header(annotated_output_fp);
  // END OF: Output annotation header

  int vcf_format_columns = 0;
  int m_dpr_column = 0;
  int m_dp_column = 0;
  int m_ro_column = 0;
  int m_dpr[100];
  int m_dp = 0;
  int m_ro = 0;
  float m_ao_percent = 0;
  char *m_pchr2;

  char *tmp;
  char m_vcf_variants[100][100];
  int m_vcf_variants_len = 0;
  char m_vcf_variant_types[100][100];
  int m_vcf_variant_types_len = 0;




  // Read variants
  while( fgets(m_temp_line, sizeof(m_temp_line), vcf_input_fp) )
  {
    if( m_temp_line[0] != '#' )
    {

      sscanf(m_temp_line, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", vcf_chr_name, &m_pos_start, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_normal_sample, vcf_sample);

      // Find all alt variants in VCF
      m_vcf_variants_len = 0;
      tmp = strtok(vcf_alt, ",");
      while( tmp != NULL )
      {
	memcpy(m_vcf_variants[m_vcf_variants_len], tmp, strlen(tmp)+1);
	m_vcf_variants_len += 1;
	tmp = strtok(NULL, ",");
      }
      // END OF: Find all alt variants in VCF


      // Find variant type from VCF info column
      tmp = strstr(vcf_info, ";TYPE=");
      if( tmp != NULL )
      {
	m_pchr2 = strtok(tmp, "=");
	m_pchr2 = strtok(NULL, ";");
	
	m_vcf_variant_types_len = 0;
	tmp = strtok(m_pchr2, ",");
	while( tmp != NULL )
	{
	  memcpy(m_vcf_variant_types[m_vcf_variant_types_len], tmp, strlen(tmp)+1);
	  m_vcf_variant_types_len += 1;
	  tmp = strtok(NULL, ",");
	}
      }
      // END OF: Find variant type from VCF info column

      // Find DPR, DP, and RO VCF format column indices
      if( m_read_variants == 0 )
      {
	m_pchr = strtok(vcf_format, ":");
	while( m_pchr != NULL )
	{
	  if( strcmp(m_pchr, g_dpr_format_str) == 0 )
	  {
	    m_dpr_column = vcf_format_columns;
	  }
	  else if( strcmp(m_pchr, g_dp_format_str) == 0 )
	  {
	    m_dp_column = vcf_format_columns;
	  }
	  else if( strcmp(m_pchr, g_ro_format_str) == 0 )
	  {
	    m_ro_column = vcf_format_columns;
	  }
	  vcf_format_columns += 1;
	  m_pchr = strtok(NULL, ":");
	}
      }
      // END OF: Find DPR, DP, and RO VCF format column indices

      // Extract DPR, DP, RO from VCF SAMPLE column
      vcf_format_columns = 0;
      m_pchr = strtok(vcf_sample, ":");
      int m_dpr_len = 0;
      char m_dpr_str[100];
      while( m_pchr != NULL )
      {
	if( vcf_format_columns == m_dpr_column )
	{
	  strcpy(m_dpr_str, m_pchr);  // Multiple values in DPR; Store DPR in string for later processing
	}
	else if( vcf_format_columns == m_dp_column )
	{
	  m_dp = atoi(m_pchr);  // Store variant sequence depth in m_dp
	}
	else if( vcf_format_columns == m_ro_column )
	{
	  m_ro = atoi(m_pchr);  // Store reference reads in m_ro
	}
	vcf_format_columns += 1;
	m_pchr = strtok(NULL, ":");
      }

      // Parse DPR values
      m_pchr = strtok(m_dpr_str, ",");
      while( m_pchr != NULL )
      {
	m_dpr[m_dpr_len] = atoi(m_pchr);
	m_dpr_len += 1;
	m_pchr = strtok(NULL, ",");
      }
      // END OF: Parse DPR values
      // END OF: Extract DPR, DP, RO from VCF SAMPLE column
	

      int m_variants_loop;
      for(m_variants_loop=0;m_variants_loop<m_vcf_variants_len;m_variants_loop++)
      {
	m_ao_percent = (float) (100*m_dpr[m_variants_loop+1]) / (float) (m_dpr[m_variants_loop+1] + m_ro);  // Percentage of reads supporting variant versus those supporting reference reads
	
	// REST/variant/variant call and find population allele frequency
	char *m_allele_freq_str;
	
	sprintf(url, "%s%s%s-%d-%s-%s", g_url_address, g_rest_variant_variant, vcf_chr_name, m_pos_start, vcf_ref, m_vcf_variants[m_variants_loop]);
	curl_easy_setopt(curl_handle, CURLOPT_URL, url);  // specify url
	curl_chunk.size = 0;  // prepare for url data
	res = curl_easy_perform(curl_handle);  // get url data
	if(res != CURLE_OK)
	{
	  fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res)); 
	  m_allele_freq_str = "-";
	}
	else  // Process REST call: Find population allele frequency
	{
	  tmp = strstr(curl_chunk.memory, "allele_freq");
	  if( tmp != NULL )
	  {
	    m_allele_freq_str = strtok(tmp, " ");
	    m_allele_freq_str = strtok(NULL, ",");
	  }
	  else
	  {
	    m_allele_freq_str = "-";
	  }
	}
	// END OF: REST/variant/variant call and find population allele frequency
	

	// REST/variant/ordered_csqs; Find most deleterious variant type/consequence
	int m_consequence_severity = m_so_terms_len;
	char m_consequence_array[1024];

	sprintf(url, "%s%s%s-%d-%s-%s", g_url_address, g_rest_variant_ordered_csqs, vcf_chr_name, m_pos_start, vcf_ref, m_vcf_variants[m_variants_loop]);
	curl_easy_setopt(curl_handle, CURLOPT_URL, url);  // specify url
	curl_chunk.size = 0;  // prepare for url data
	res = curl_easy_perform(curl_handle);  // get url data
	if(res != CURLE_OK) 
	{
	  fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
	  memcpy(m_consequence_array, "-", 2 );
	}
	else
	{
	  for( m_a_loop=0;m_a_loop<m_so_terms_len;m_a_loop++)
	  {
	    if( strstr(curl_chunk.memory, m_so_terms[m_a_loop]) != NULL )
	    {
	      m_consequence_severity = m_a_loop;
	      break;
	    }
	  }
	  if( m_consequence_severity >= m_so_terms_len )
	  {
	    memcpy(m_consequence_array, "-", 2 );
	  }
	  else
	  {
	    memcpy(m_consequence_array, m_so_terms_no_quotes[m_consequence_severity], strlen(m_so_terms_no_quotes[m_consequence_severity])+1);
	  }
	}
	// END OF: REST/variant/ordered_csqs; Find most deleterious variant type/consequence
	
	// Output annotation
	fprintf(annotated_output_fp, "%s\t%d\t%s\t%s\t%d\t%d\t%f\t%s\t%s\t%s\n", vcf_chr_name, m_pos_start, vcf_ref, m_vcf_variants[m_variants_loop], m_dp, m_dpr[m_variants_loop+1], m_ao_percent, m_allele_freq_str, m_vcf_variant_types[m_variants_loop], m_consequence_array );
	// END OF: Output annotation
	
	m_read_variants += 1;
	if( m_read_variants / 100 * 100 == m_read_variants || m_read_variants == 1 )
	{
	  printf("\nVariants processed: %d", m_read_variants);
	}
      }
    }
  }
  printf("\nRead %d variants from %s\n\n", m_read_variants, vcf_input_file_name);
  // END OF: Read variants


  // Curl cleanup
  curl_easy_cleanup(curl_handle);
  curl_global_cleanup();
  // END OF: Curl cleanup


  // Free memory
  free(curl_chunk.memory);
  free_2d_str_array(m_so_terms_no_quotes, m_so_terms_len);
  free_2d_str_array(m_so_terms, m_so_terms_len);
  // END OF: Free memory


  // Close files
  fclose(vcf_input_fp);
  fclose(annotated_output_fp);
  // END OF: Close files

  return 0;  

}


