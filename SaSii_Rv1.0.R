#-------REQUESTS-------------------------

#----------------------------------------

#-------LIBRARIES------------------------

#----------------------------------------

# ------SOME GLOBAL VARIABLES------------
args <- commandArgs(trailingOnly = TRUE);

os = Sys.info()['sysname'];
path = '';
lines_array <- vector(); #VERIFICAR ESSA VARIAVEL
line = '';
input_lines <- vector();
local_multi_he = 0; # temp variable for storage multilocus He.
local_multi_ho = 0; # temp variable for storage multilocus Ho.

localtime = format(Sys.time(), "%c");

less_frq_all = more_frq_all = 0;
less_frq_allname = more_frq_allname = '';
freqMin = 0;
he_diff = ho_diff = list();
n_size = 0;
#----------------------------------------

#----SUBROTINES--------------------------
print_hashes <- function(target_hash1){
    sorted_target_hash = sort(names(target_hash1));
    for(locus in sorted_target_hash){
        print(paste("Locus =", locus));
        sorted_locus = sort(target_hash1[[locus]]);

        for(key in sorted_locus)
            print(paste(key, "=", sprintf("%13.10f", target_hash[[locus]][key])));
    }
}
#----------------------------------------
print_HoA <- function(HoA){
    sorted_HoA = sort(names(HoA));

    for(locus in sorted_HoA){
        print(paste("Locus =", locus));
        print(paste(HoA[[locus]])); 
    }
}
#----------------------------------------
print_HoHoA <- function(HoHoA){
    if (HoHoA == ''){
        print("HoHoA not defined !");
        quit();
    }
    print(paste("HoHoA:", HoHoA));

    sorted_HoHoA = sort(names(HoHoA));

    for(local_n in sorted_HoHoA){
        print(paste("N =", local_n));
        sorted_local_n = sort(names(local_n));

        for(array in sorted_local_n){
            print(paste("array =", array));
            print(HoHoA[[local_n]][array]);
        }
    }
}
#----------------------------------------
save_array_output <- function(lines, mid_name){

    output_name = paste0(mid_name,'.sso');
    #output_name =~ s/.csv/_/g; VERIFICAR COMO FAZER ISSO

    path = "results/";
    file.create(paste0(path, output_name));

    fileConn<-file(paste0(path, output_name));
    writeLines(lines, fileConn);
    close(fileConn);
}
#----------------------------------------
save_hash_output <- function(lines){
    #FALTA ENTENDER ESSA FUNÇÃO
}
#----------------------------------------
print_hashes_to_arrays <- function(locus, target_hash1){
    #FALTA ENTENDER ESSA FUNÇÃO
}
#----------------------------------------
# Handle csv files
csv_adjust <- function(csv_data){ #VERIFICAR COMO PODEM VIR OS CSVs PARA SABER COMO TRATAR OS DADOS
    locus = "";
    csv_data = csv_data[-length(csv_data)];
    for(lines in names(csv_data)){
        if(grepl('X', lines, fixed = TRUE)) {
            names(csv_data)[names(csv_data) == lines] <- paste0(locus, '.2'); #verificar isso depois
        }else{locus = lines;}
    }
    return (csv_data);
}
#----------------------------------------
handleAllele_Gnot <- function(auxAlleles, allele){
  if(allele %in% names(auxAlleles)){
    auxAlleles[allele] = as.integer(auxAlleles[allele]) + 1;
  }else{
    if(is.na(auxAlleles)){
      auxAlleles = array(1);
    }else{
      auxAlleles <- append(auxAlleles, 1);
    }
    names(auxAlleles)[length(auxAlleles)] = allele;
  }
  return(auxAlleles);
}
#----------------------------------------
make_hash_count <- function(data){
    #0,18 seg
    # alleles_unsorted = vector(); # Create array to temp storage of unsorted alleles
    # alleles_sorted = vector();	# Create array to temp storage of sorted alleles
    # #print("1");
    # for(indv in 1:nrow(data)){ #For each individual
    #     for(pos in seq(1, length(data[indv,]), ploid)){ #Select the first position (column) of each locus in data.frame
    #         auxAlleles = allele_numbers[[loci_names[(pos+1)/ploid]]]; #Get current state of locus
    #         #print("2");
    #         #Count alleles in locus for each individual
    #         for(allele_pos in pos:(pos+ploid-1)){
    #             index = as.character(data[indv, allele_pos]);
    #             #print("3");
    #             if (data[indv, allele_pos] %in% names(auxAlleles)){
    #                 auxAlleles[index] = as.integer(auxAlleles[index]) + 1;
    #             }else{
    #                 if(is.na(auxAlleles[1])){ #Verify if array is empty and inicialize it
    #                     auxAlleles = array(1); 
    #                     names(auxAlleles) = index;
    #                 }
    #                 else{ #Concatenates a new allele for counting
    #                     auxAlleles = append(auxAlleles, 1);
    #                     names(auxAlleles)[length(auxAlleles)] = index;
    #                 }
    #             }
    #         }
    #         #print("4");
    #         #print(typeof(auxAlleles));
    #         allele_numbers[loci_names[(pos+1)/ploid]] <<- list(auxAlleles); #Updates the allele lists

    #         #Generate and count each genotype
    #         #if(teste) print(sort(unlist(data[indv, pos:(pos+ploid-1)])));
    #         #if(teste) quit();
    #         auxGenotypes = gnot_numbers[[loci_names[(pos+1)/ploid]]]; #Get current state of genotypes
    #         allele_locus = sort(unlist(data[indv, pos:(pos+ploid-1)])); #Get all alleles(sorted) on this locus from this individual
    #         #print("5");
    #         #if(teste) quit();
    #         for(allele1 in allele_locus){
    #             allele_locus = allele_locus[-1];
    #             for(allele2 in allele_locus){
    #                 genotype = paste0(allele1, "-", allele2);
    #                 #print("6");
    #                 if (genotype %in% names(auxGenotypes)){
    #                     auxGenotypes[genotype] = as.integer(auxGenotypes[genotype]) + 1;
    #                 }else{
    #                     if(is.na(auxGenotypes[1])){ #Verify if array is empty and inicialize it
    #                     auxGenotypes = array(1); 
    #                     names(auxGenotypes) = genotype;
    #                     }
    #                     else{ #Concatenates a new allele for counting
    #                     auxGenotypes = append(auxGenotypes, 1);
    #                     names(auxGenotypes)[length(auxGenotypes)] = genotype;
    #                     }
    #                 }
    #             }
    #         }
    #         gnot_numbers[loci_names[(pos+1)/ploid]] <<- list(auxGenotypes); #Updates the allele lists
    #     }
    # }
    #start = Sys.time();
    #0.09
    # a = lapply(1:nrow(data), function(line){
    #     data_line = data[line,];
    #     lapply(loci_names, function(locus){
    #         allele1 = as.character(data_line[1]);
    #         allele2 = as.character(data_line[2]);
    #         data_line <<- data_line[-c(1,2)];

    #         auxAlleles = handleAllele_Gnot(allele_numbers[[locus]][[1]], allele1);
    #         auxAlleles = handleAllele_Gnot(auxAlleles, allele2);
    #         allele_numbers[[locus]] <<- list(auxAlleles);

    #         sorted_allele = sort(c(allele1, allele2));
    #         genotype = paste0(sorted_allele[1], "-", sorted_allele[2]);
    #         auxGnot = handleAllele_Gnot(gnot_numbers[[locus]][[1]], genotype);
    #         gnot_numbers[[locus]] <<- list(auxGnot);
    #     });
    # });

    #0.03
    a = lapply(loci_names, function(locus){
        alleles = unlist(subset(data, TRUE, c(locus, paste0(locus,".2"))));
        #print(alleles);
        counted_alleles = table(alleles);
        allele_numbers[[locus]] <<- list(counted_alleles); 

        counted_gnot = table(paste0(alleles[seq(1, length(alleles)/2)], "-", alleles[seq(length(alleles)/2 + 1, length(alleles))]))
        gnot_numbers[[locus]] <<- list(counted_gnot);
    });
    #print(allele_numbers)
    #0,13 seg
    # for (line in 1:nrow(data)) {
    #     data_line = data[line,];
    #     for(locus in loci_names){
    #         allele1 = as.character(data_line[1]);
    #         allele2 = as.character(data_line[2]);
    #         data_line = data_line[-c(1,2)];

    #         auxAlleles = allele_numbers[[locus]][[1]];
    #         auxAlleles = handleAllele_Gnot(auxAlleles, allele1);
    #         auxAlleles = handleAllele_Gnot(auxAlleles, allele2);
    #         allele_numbers[[locus]] <<- list(auxAlleles);

    #         sorted_allele = sort(c(allele1, allele2));
    #         genotype = paste0(sorted_allele[1], "-", sorted_allele[2]);
    #         auxGnot = handleAllele_Gnot(gnot_numbers[[locus]][[1]], genotype);
    #         gnot_numbers[[locus]] <<- list(auxGnot);
    #     }
    # }
    #end = Sys.time();
    #print(paste0("Tempo make_hash: ", (end-start)));
    #quit();
    #print("");
}
#----------------------------------------
allele_count <- function(allele_hash){
    count = 0;
    a = lapply(allele_hash, function(locus){
        lapply(names(locus[[1]]), function(allele){
            if(allele != '0'){
                count <<- count + 1;
            }
        })
    });
    # for(locus in allele_hash){
    #     for(allele in names(locus[[1]])){
            
    #     }
    # }
    return(count);
}
#----------------------------------------
allelic_richness <- function(allele_hash){ #VERIFICAR ESSA FUNÇÃO
    a = lapply(names(allele_hash), function(locus){
        allele_rich[[locus]] <<- length(allele_hash[[locus]]); #em cada locus precisa armazenar a quantidade de alelos existentem no hash
    });
}
#----------------------------------------
freq_calc <- function(numbers_hash){
    freq_hash = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(freq_hash) <- loci_names;
    sorted_numbers_hash = sort(names(numbers_hash));

    a = lapply(sorted_numbers_hash, function(locus){
        sorted_locus_hash = sort(names(numbers_hash[[locus]][[1]]));
        #Alleles is a list without lost data
        good_data = as.character(unlist(
            lapply(sorted_locus_hash, 
                   function(x) x[(which(x != '000' & x != '000-000' & 
                                        x != '0' & x != '0-0'))])));

        freq_hash[[locus]] <<- list(numbers_hash[[locus]][[1]][good_data]); #Creates a copy of the variable type
        total_alleles = sum(numbers_hash[[locus]][[1]][good_data]);

        freq_hash[[locus]] <<- numbers_hash[[locus]][[1]][good_data] / total_alleles; #Calculate the frequency for each element
    });
    return(freq_hash);
}
#----------------------------------------
find_less_more <- function(freqs){
    #Extract the position of the less freq allele
    less_pos = match(min(unlist(freqs)), unlist(freqs))
    less_locus = strtrim(names(unlist(freqs)[less_pos]), 6);
    less_allele = substr(names(unlist(freqs)[less_pos]), 8, nchar(names(unlist(freqs))));

    #Extract the position of the more freq allele
    more_pos = match(max(unlist(freqs)), unlist(freqs))
    more_locus = strtrim(names(unlist(freqs)[more_pos]), 6);
    more_allele = substr(names(unlist(freqs)[more_pos]), 8, nchar(names(unlist(freqs))));

    local_less = paste0(less_locus, "-", less_allele);
    local_more = paste0(more_locus, "-", more_allele);
    return (list(local_less, local_more));
}
#----------------------------------------
less_more_data <- function(freqs){
    aux_array = strsplit(less_frq_all, '-')[[1]];
    locus_l = aux_array[1]; allele_l = aux_array[2];
    aux_array = strsplit(more_frq_all, '-')[[1]];
    locus_m = aux_array[1]; allele_m = aux_array[2];

    local_meanL = local_SDL = local_meanM = local_SDM = 0;
    less_min = less_max = more_min = more_max = 0;

    #Sort the frequencies of the locus_alleles
    less_locus_freqs = sort(freqs[[locus_l]][[allele_l]][-1]);
    more_locus_freqs = sort(freqs[[locus_m]][[allele_m]][-1]);

    #Get frequencies of less and more locus_alleles
    less_min = less_locus_freqs[1];
    less_max = less_locus_freqs[2];
    more_min = more_locus_freqs[1];
    more_max = more_locus_freqs[2];

    #Calculate the mean and SD of locus
    local_meanL = sprintf("%.5f", mean(less_locus_freqs));
    local_meanM = sprintf("%.5f", mean(more_locus_freqs));
    local_SDL = sprintf("%.5f", sd(less_locus_freqs));
    local_SDM = sprintf("%.5f", sd(more_locus_freqs));

    # sorted_freqs = sort(names(freqs));

    # for(locus in sorted_freqs){
    #     sorted_locus = sort(names(freqs[[locus]]));
    
    #     for(allele in sorted_locus){
    #         arrays_freqs = sort(names(freqs[[locus]][allele]));

    #         #VERIFICAR COMO USAR SHIFT E POP EM R
    #         if(locus == locus2 && allele == allele2){
    #             aux_array = mean_and_sd(arrays_freqs);
    #             local_meanL = aux_array[1]; local_SDL = aux_array[2];

    #             less_min = shift(arrays_freqs);
    #             less_max = pop(arrays_freqs);
    #         }

    #         if(locus == locus3 && allele == allele3){
    #             aux_array = mean_and_sd(arrays_freqs);
    #             local_meanM = aux_array[1]; local_SDM = aux_array[2];

    #             more_min = shift(arrays_freqs);
    #             more_max = pop(arrays_freqs);
    #         }
    #     }
    # }

    local_data = vector();

    less_min = sprintf('%.5f', less_min);
    less_max = sprintf('%.5f', less_max); 
    local_data <- append(local_data, paste0(n_size, "\t", paste0("lessfreq", less_frq_all), "\t", local_meanL, "\t", local_SDL, "\t", less_min, "\t", less_max));

    more_min = sprintf('%.5f', more_min);
    more_max = sprintf('%.5f', more_max);
    local_data <- append(local_data, paste0(n_size, "\t", paste0("morefreq", more_frq_all), "\t", local_meanM, "\t", local_SDM, "\t", more_min, "\t", more_max));

    return(local_data);
}
#----------------------------------------
expct_GNOT_freq <- function(alelic_freq){ # Calculate expected genotypes frquencies from HoH with allelic freq
    EXPCT_GNOT_freq = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(EXPCT_GNOT_freq) <- loci_names; # Hash with expected genotypes frequencies

    sorted_alelic_freq = sort(names(alelic_freq));

    a = lapply(sorted_alelic_freq, function(locus){
        exp_gnot_freq = 0;
        generated_genotypes = generate_genotypes(names(alelic_freq[[locus]]));

        EXPCT_GNOT_freq[locus] = list(array(rep(array(), length(generated_genotypes)), c(length(generated_genotypes), 1)));
        names(EXPCT_GNOT_freq[[locus]]) <- generated_genotypes;

        lapply(generated_genotypes, function(genotype){
            
            aux_array = strsplit(genotype, '-')[[1]]; #VERIFICAR ESSA PARTE
            allele1 = aux_array[1]; allele2 = aux_array[2];

            if(allele1 == allele2){
                exp_gnot_freq = alelic_freq[[locus]][allele1] * alelic_freq[[locus]][allele2];
            }else {
                if (allele1 != allele2) {
                    exp_gnot_freq = 2 * (alelic_freq[[locus]][allele1] * alelic_freq[[locus]][allele2]);
                }else { #Acho que isso nunca aconteceria
                    exp_gnot_freq = "Wrong !!!";
                    print(paste(allele1, "fail to compares to", allele2));
                }
            }

            EXPCT_GNOT_freq[[locus]][genotype] = exp_gnot_freq;
        });
    });

    return(EXPCT_GNOT_freq);
}
#----------------------------------------
generate_genotypes <- function(allelic_freq_inlocus) { # Create list of possible genotypes from an hash: allele => value
    generated_genotypes = array(); # Array to save genotypes list

    sorted_allelic_freq_inlocus = sort(allelic_freq_inlocus); #VERIFICAR

    lapply(sorted_allelic_freq_inlocus, function(allele1){
        lapply(sorted_allelic_freq_inlocus, function(allele2) {
            genotype = paste0(allele1, '-', allele2);

            if(!is.na(generated_genotypes[1])){
                generated_genotypes <<- append(generated_genotypes, genotype);
            }else{
                generated_genotypes <<- array(genotype);
            }
        });
        sorted_allelic_freq_inlocus <<- sorted_allelic_freq_inlocus[-1]; #remove o primeiro elemento da lista
    });
    return(generated_genotypes);
}
#----------------------------------------
ho_calc <- function(gnot_freq) { # Find observed heterosygosity for data frame
    sorted_gnot_freq = sort(names(gnot_freq));

    a = lapply(sorted_gnot_freq, function(locus){
        sorted_locus = sort(names(gnot_freq[[locus]]));
        ho[[locus]] <<- 0;
        
        lapply(sorted_locus, function(genotype){
            aux_array = strsplit(genotype, '-')[[1]];
            allele1 = aux_array[1]; allele2 = aux_array[2];
            
            if(allele1 != allele2)
                ho[[locus]] <<- ho[[locus]] + gnot_freq[[locus]][genotype];
        });
    });
}
#----------------------------------------
ho_range <- function(local_ho) {
    sorted_local_ho = sort(names(local_ho));
    a = lapply(sorted_local_ho, function(locus){
        ho_for_range[[locus]] <<- append(ho_for_range[[locus]], ho[[locus]]);
    });#ho_for_range = ho; #VERIFICAR ISSO
    sorted_ho = sort(names(ho));

    ho_array_local = vector();
    a = lapply(sorted_ho, function(locus){
        ho_array_local <<- append(ho_array_local, ho[[locus]]);
    });
    
    #aux_array = mean_and_sd(ho_array_local);
    mean = sprintf("%.5f", mean(ho_array_local)); SD = sprintf("%.5f", sd(ho_array_local));

    local_multi_ho <<- mean;
}
#----------------------------------------
random_sampling <- function(n) {
    count = n; 
    total_lines = total_individuals;
    sample = vector();
    backup_input <- csv_data;
    #teste = 2;
    while(count > 0) {
        #print(paste("Tamanho:",nrow(backup_input)));
        number = as.integer(runif(1, min=1, max=nrow(backup_input)));
        #print(number);
        #Creates a matrix with the random individuals selected
        if(count == n){
            sample = matrix(data = backup_input[number,], nrow = 1, ncol = 18);
        }else{
            #sample = rbind(sample, as.list(backup_input[-2,]));
            sample = rbind(sample, as.list(backup_input[number,]));
        }
        #sample = matrix(data = data[1,], nrow=1, ncol=18)
        #sample = rbind(sample, as.list(data[4,]))
        #print(sample);
        #push(sample, backup_input[number,]); #Extracts randomized line from matrix
        #length = total_lines - number - 1;

        backup_input = backup_input[-number,]; #Removes the specific line from matrix

        #if (length > 0){ backup_input = backup_input[-c(number:number+length),];}
        #else{ backup_input = backup_input[-1,];}

        #total_lines = total_lines - 1;
        count = count - 1;
        #teste = teste + 1;
    }
    #print(sample);
    #quit();
    
    return(sample);
}
#----------------------------------------
sample_calcs <- function(sample) {
    allele_rich <<- gnot_numbers <<- allele_numbers <<- array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(allele_rich) <<- names(gnot_numbers) <<- names(allele_numbers) <<- loci_names; #cria o nome da coluna

    make_hash_count(sample);
    allelic_richness(allele_numbers);
    allele_freq <<- freq_calc(allele_numbers);
    gnot_freq = freq_calc(gnot_numbers);
    EXPCT_GNOT_freq = expct_GNOT_freq(allele_freq);
    he_calc(allele_freq, n_size, '2');
    multi_he = local_multi_he;

    he_mean <<- append(he_mean, as.numeric(multi_he));
    he_diff <<- het_diff_calc(he_orig, he_diff, he);

    mean_freq_prep1(allele_freq);

    sorted_orig_gnot_freq = sort(names(orig_gnot_freq));

    #for (locus in sorted_orig_gnot_freq) # Inicialize %Ho with all loci
    #    ho[[locus]] = 0; #VERIFICAR SE ISSO ESTA OK
    #print(paste0("Class: ", n_size, "    Repeat N: ", repeats));
    ho_calc(gnot_freq);
    ho_range(ho);
    
    #print(local_multi_ho);
    
    multi_ho = local_multi_ho;
    ho_mean <<- append(ho_mean, as.numeric(multi_ho));

    ho_diff <<- het_diff_calc(ho_orig, ho_diff, ho);
    nei_fst <<- calc_nei_fst(multi_ht, multi_he_orig, multi_he, n_size);

    local_nei_d = nei_d(allele_freq);
    neiD_array <<- local_nei_d;

    local_rogers_d = rogers_d(allele_freq);
    rogersD_array <<- local_rogers_d;
}
#----------------------------------------
mixed_calcs <- function(sample) {
    allele_rich <<- gnot_numbers <<- allele_numbers <<- array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(allele_rich) <<- names(gnot_numbers) <<- names(allele_numbers) <<- loci_names; #cria o nome da coluna

    make_hash_count(sample);
    allele_freq <<- freq_calc(allele_numbers);
    
    gnot_freq = freq_calc(gnot_numbers);
    EXPCT_GNOT_freq = expct_GNOT_freq(allele_freq);
    he_calc(allele_freq, n_size, '1');

    mix_he = he;
    #he = ();
    multi_ht <<- local_multi_he;
    return(mix_he);
}
#----------------------------------------
he_calc <- function(local_hash, number_id, type) { #Needs a HoH with allelic frequencies
    somatory = 0;
    r = 0; #Number o loci

    sorted_local_hash = sort(names(local_hash));

    a = lapply(sorted_local_hash, function(locus){
        somatory = sum(local_hash[[locus]] ** 2);

        he_simple = 1 - somatory;
        he[[locus]] <<- he_simple * ((2*number_id)/((2*number_id)-1));
    
        if(type == '2') {
            # print(typeof(he_for_range[[locus]]));
            # print(he[[locus]]);
            # print(append(he_for_range[[locus]], he[[locus]]));
            he_for_range[[locus]] <<- append(he_for_range[[locus]], he[[locus]]); 
            #push(he_for_range[[locus]], he[[locus]]); #VERIFICAR ISSO
        }
    });
    #print(he_for_range);

    local_multi_he <<- 0;
    he_array_local = vector();
    sorted_he = sort(names(he));

    a = lapply(sorted_he, function(locus){
        he_array_local <<- append(he_array_local, he[[locus]]);
    });
    #aux_array = mean_and_sd(he_array_local);

    mean = sprintf("%.5f", mean(he_array_local)); SD = sprintf("%.5f", sd(he_array_local));
    local_multi_he <<- mean;
}
#----------------------------------------
calc_nei_fst <- function(ht, hsA, hsB, nB) { #Calcs Nei Fst (1973). Needs $Multi_Ht, $Multi_He_ORIGINAL, $Multi_He, $n_size.
    nA = total_individuals;
    nei_fst = 0;
    
    upper1a = ((nA * as.numeric(hsA)) + (as.numeric(nB) * as.numeric(hsB))) / (nA + as.numeric(nB));
    upper1b = as.numeric(ht) - upper1a;
    nei_fst = upper1b / as.numeric(ht);

    return (nei_fst);
}
#----------------------------------------
#TALVEZ NÃO PRECISE USAR ISSO
sample_for_R <-function(sample)  #clean id names for use with Geneland package in R
{
  new_sample = itens_array = vector();
  
  for (line in 1:nrow(sample)){
    itens_array = strsplit(sample[line,], ' ')[[1]];
    itens_array = itens_array[-1]; 
    for (item in itens_array) #verificar se isso funciona igual ao perl
        line = paste0(line, ' ', item); #VERIFICAR ISSO

    new_sample = append(new_sample, line);
  }
  return ( new_sample );
}
#----------------------------------------
nei_d <- function(freq2) { #Nei distance (1972)
    freq1 = orig_allele_freq;

    jx = jy = jxy = l = 0; 
    jx_array = jy_array = jxy_array = vector();

    sorted_freq1 = sort(names(freq1));

    lapply(sorted_freq1, function(locus) {
        l = l + 1;

        sorted_locus = sort(names(freq1[[locus]]));

        lapply(sorted_locus, function(allele) {
            jx = freq1[[locus]][allele] ** 2;
            jx_array <- append(jx_array, jx);

            if ( is.na(freq2[[locus]][allele]) ){
                freq2[[locus]][allele] = 0;
            }
            jy = freq2[[locus]][allele] ** 2;
            jy_array <- append(jy_array, jy);
        });
    });

    lapply(sorted_freq1, function(locus) {
        sorted_locus = sort(names(freq1[[locus]]));

        lapply(sorted_locus, function(allele) {
            if(is.na(freq2[[locus]][allele])){ jxy = 0;}
            else {
                jxy = freq1[[locus]][allele] * freq2[[locus]][allele];
            }
            jxy_array <- append(jxy_array, jxy);
            #print(paste0(freq1[[locus]][allele], " * ",freq2[[locus]][allele], " = ",jxy));
        });
    });

    jx = sum(jx_array)/l;
    jy = sum(jy_array)/l;
    jxy = sum(jxy_array)/l;

    gx = jx;
    gy = jy;

    neiD = jxy / sqrt(gx * gy);
    neiD = log(neiD);
    neiD = 0 - neiD;

    return(neiD);
}
#----------------------------------------
rogers_d <- function(freq2) { # Calculares modified Rogers distance. Needs HoH of resample allele frequencies. Roger (1972), Wright (1978), GOodman and Stuber (1983) 

    freq1 = orig_allele_freq;

    diff_sum = 0; l = 0;
    diff_sum_array = vector();

    sorted_freq1 = sort(names(freq1));

    lapply(sorted_freq1, function(locus){
        l = l + 1;
        diff_array = vector();
        #print(locus);
        sorted_locus = sort(names(freq1[[locus]]));

        lapply(sorted_locus, function(allele) {
            if(is.na(freq2[[locus]][allele])) freq2[[locus]][allele] = 0; #Treat cases where allele do not exist
            
            diff = freq1[[locus]][allele] - freq2[[locus]][allele];
            diff = diff ** 2;
            
            diff_array <- append(diff_array, diff);
        });

        diff_sum = sum(diff_array);
        diff_sum_array <- append(diff_sum_array, diff_sum);
    });
    #print(diff_sum_array);
    rogersD = sqrt(sum(diff_sum_array));
    rogersD = rogersD / sqrt(2 * l);

    return (rogersD);
}
#----------------------------------------
mean_freq_prep1 <- function(freqs) {
    sorted_orig_allele_freq = sort(names(orig_allele_freq));

    for(locus in sorted_orig_allele_freq){
        sorted_locus = sort(names(orig_allele_freq[[locus]]));

        for(allele in sorted_locus){
            #print(paste(locus," - ",allele, "= ", allele_freq[[locus]][allele]));
            if(is.na(allele_freq[[locus]][allele])){
                allele_freq[[locus]][allele] <<- 0;
                #print(allele_freq[[locus]][allele]);
                #print(freqs[[locus]][allele]);
            }
            mean_allele_freq[[locus]][allele] <<- list(append(unlist(mean_allele_freq[[locus]][allele]), allele_freq[[locus]][allele]));
        }
    }
    #print(mean_allele_freq);
    #quit();
}
#----------------------------------------
freq_diff_calc <- function(freqs) {
    mean_freqs = freqs;
    sorted_freqs = sort(names(freqs));
    freqs_differences = vector();

    for(locus in sorted_freqs) {
        sorted_locus = sort(names(freqs[[locus]]));

        for(allele in sorted_locus) { #VERIFICAR SE USAR ESSA FUNÇÃO NÃO DARIA PROBLEMAS
            mean_freqs[[locus]][allele] = sprintf('%.5f', sum(freqs[[locus]][[allele]][-1])/length(freqs[[locus]][[allele]][-1]));

            diff = orig_allele_freq[[locus]][allele] - as.double(mean_freqs[[locus]][[allele]]);
            diff = abs(diff);

            if(allele < 100) {allele = paste0('_', allele);}

            freqs_differences <- append(freqs_differences, paste0(n_size, "\t", paste0(locus, ".", allele), "\t", sprintf('%.5f', diff)));
        }
    }

    #print (freqs_differences[2]);
    #quit();
    #PARTE DO CÓDIGO JUNTADO 
    # sorted_orig_allele_freq = sort(names(orig_allele_freq));

    # for(locus in sorted_orig_allele_freq){
    #     sorted_locus = sort(orig_allele_freq[[locus]]);

    #     for(allele in sorted_locus){
    #         print(orig_allele_freq[[locus]][allele]);
    #         print(mean_freqs[[locus]][allele]);
    #         quit();
    #         #diff = orig_allele_freq[[locus]][allele] - mean_freqs[[locus]][allele];
    #         diff = abs(diff);

    #         if(allele < 100) {allele = paste0('_', allele);}

    #         diff = sprintf('%.5f', diff); #VERIFICAR ISSO
    #         freqs_differences <- append(freqs_differences, paste(n_size, "\t", paste0(locus, allele), "\t", diff, "\n"));
    #     }
    #     print (freqs_differences);
    #     quit();
    # }
    #return ();
    return (freqs_differences);
}
#----------------------------------------
freq_5 <- function(freq_hash) {
    five_percent_hash = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(five_percent_hash) <- loci_names; # HoH to save alleles with freq higher than freqMin
    five_percent_n = 0;

    sorted_freq_hash = sort(names(freq_hash));

    for(locus in sorted_freq_hash) {
        sorted_locus = sort(names(freq_hash[[locus]]));

        for(key in sorted_locus){
            if(as.numeric(freq_hash[[locus]][key]) >= freqMin){
                if(is.na(five_percent_hash[[locus]])){
                    five_percent_hash[[locus]] = c(key);
                }else{
                    five_percent_hash[[locus]] = list(c(unlist(five_percent_hash[[locus]]), key));
                }
                five_percent_n = five_percent_n + 1; #writes in local variable with same name as global one
            }
        }
    }

    return(list(five_percent_hash, five_percent_n)); #returning hash as reference 
}
#----------------------------------------
count_freq_5_sample <- function(freq_hash, sample_local) {
    reference_hash = orig_five_percent;
    local_ref_count = five_percent_n;
    general_local_count = 0;

    sorted_freq_hash = sort(names(freq_hash));
    
    for(locus in sorted_freq_hash) {
        local_array = as.array(unlist(reference_hash[[locus]]));
        sorted_freq_locus = sort(names(freq_hash[[locus]]));
        local_count = 0;

        for(key in sorted_freq_locus){
            #print(paste(key, "= ", freq_hash[[locus]][key], key %in% local_array));
            if((freq_hash[[locus]][key] != 0) && (key %in% local_array)) { #verificar se está correto
                local_count = local_count + 1;
                general_local_count = general_local_count + 1;
            }
        }
        local_count = local_count / length(local_array);
        sample_local[[locus]] <- append(sample_local[[locus]], local_count);
    }
    #print(paste(general_local_count,"/",local_ref_count));
    general_local_count = general_local_count / local_ref_count;

    return (list(general_local_count, sample_local));
}
#----------------------------------------
r_input_HoH_HoA <- function(target_hash) {
    sorted_target_hash = sort(names(target_hash));
 
    for(locus in sorted_target_hash) {
        sorted_locus = sort(names(target_hash[[locus]]));

        for(key in sorted_locus) {
            array_name = paste0(locus, '_', key);

            if(is.null(freq_arrays_complete[[as.character(n_size)]][[array_name]])){
                freq_arrays_complete[[as.character(n_size)]][[array_name]] <<- list(target_hash[[locus]][key]);
            }else{
                freq_arrays_complete[[as.character(n_size)]][[array_name]] <<- list(append(unlist(freq_arrays_complete[[as.character(n_size)]][[array_name]]), target_hash[[locus]][key]));
            }
        }
    }
}
#----------------------------------------
r_input_HoHoA_array <- function(target_hash1) { #VERIFICAR POIS AQUI ELE DIZ QUE USA COMO REFERENCIA
    array_final = list(paste0('N', '\t', 'ID', '\t', 'repeat_values'));
    sorted_target_hash1 = sort(names(target_hash1));
    print(names(target_hash1));
    for(local_n in sorted_target_hash1) {
        sorted_local_n = sort(names(local_n));
        print(sorted_local_n);
        for(array_ in sorted_local_n) {
            print(array_final);
            array_final <- append(array_final, paste0(local_n, '\t', array_, '\t', target_hash1[[local_n]][[array_]]));
        }
    }
    print(array_final);
    return(array_final);
}
#----------------------------------------
make_path_files <- function(localfile, n){
    file.create(paste0("fig", n, ".dat"));
    #file.create(paste0("results/", localfile,".sso"));
    dput(paste0("results/", localfile,".sso"), 
         file = paste0("fig", n, ".dat"));
}
#----------------------------------------
make_temp_files <- function() { #PRECISA FINALIZAR AINDA
    n = 0;
    for(n in c(1:1)){
        file_name = paste0("fig", n, ".dat");
        out_file = "";
        i = j = q = w = N = garb = 0;
        ind_n = locus_n = linhas = n_linha = 1;
        matriz = vector(); aux = "?"; str_aux = "";

        file_data <- dget(file = file_name);

        print(paste("Arquivo aberto para leitura:", file_name));

        file_name = paste0(file_name, ".out");

    }
}
#----------------------------------------
plot_graphs <- function() {
    source("SaSii_plots.R");
}
#----------------------------------------
het_diff_calc <- function(hash, local_hash, orig_hash) {
    sorted_orig_hash = sort(names(orig_hash));
    #print("HET DIFF");
    #print(abs(orig_hash[["MpeB03"]] - hash[["MpeB03"]]));
    for(locus in sorted_orig_hash) {
        #print(orig_hash[[locus]]);
        #print(local_hash[[locus]]);
        #print(paste(orig_hash[[locus]], "-", hash[[locus]], "=", abs(orig_hash[[locus]] - hash[[locus]])));
        #print(append(local_hash[[locus]], (abs(orig_hash[[locus]] - hash[[locus]]))));
        local_hash[[locus]] <- append(local_hash[[locus]], abs(orig_hash[[locus]] - hash[[locus]]));
        #print(local_hash[[locus]]);
    }
    # print(orig_hash[["MpeB03"]]);
    # print(hash[["MpeB03"]]);
    # print(abs(orig_hash[["MpeB03"]] - hash[["MpeB03"]]));
    # print(append(local_hash[["MpeB03"]], (abs(orig_hash[["MpeB03"]] - hash[["MpeB03"]]))));
    # print("Total");
    # print(local_hash[["MpeB03"]]);
    return(local_hash);
}
#----------------------------------------
het_diff_push <- function(hash) {
    array = vector();
    sorted_hash = sort(names(hash));

    for(locus in sorted_hash) {
        hash[[locus]] <- hash[[locus]][-1]; # take of '5' from beggining of hash VERIFICAR ISSO

        het_men_diff = sprintf("%.5f", mean(hash[[locus]])); het_sd_diff = sd(hash[[locus]]);
        array <- append(array, paste0(n_size, "\t", locus, "\t", het_men_diff, "\t", ifelse(is.na(het_sd_diff), 0, sprintf("%.5f", het_sd_diff))));
    }

    return (array);
}
#----------------------------------------

#----MAIN CODE---------------------------
print("SaSii - Sample_Size Impact , 'R Beta version ', (2021)");
print(paste("Running in:",os));

if(length(args) == 0){
    cat("Put you input file path: ");
    path = as.character(readLines("stdin", n=1));
}else{
    path = args;
}

if(path == "help" | path == "ajuda"){
    print("Help file not exist at this version.\nDone.");
    quit();
} else{
    if(path == "exit") {
        print("Done."); quit();
    } else {
        if(path == '' && path == '\n'){
            print("WARNING: You passed no input file.\n Try: sample-size-simulationsXX.pl inputfile.txt");
            quit();
        }
    }
}



#Check if results folder exists, if do, delete it.
if(dir.exists("results")){
    unlink("results", recursive = TRUE);
}
dir.create("results");

#Check if exists any temp file and deletes it.
if(!is.null(list.files(paste(pattern = ".dat")))){
    print("Cleaning temp files....");
    unlink(paste0("fig","*",".dat"));
}

make_path_files('1-5percent_rate', 1);
make_path_files('2-freq_dif', 2);
make_path_files('3-freq_impact', 3);
make_path_files('4-He_impact', 4);
make_path_files('5-meanHe_impact', 5);
make_path_files('6-Fst', 6);
make_path_files('7-Nei', 7);
make_path_files('8-Roger', 8);
make_path_files('9-meanHo_impact', 9);
make_path_files('10-Ho_impact', 10);
make_path_files('11-Ho_diff', 11);
make_path_files('12-He_diff', 12);


# Defining the begining of data.---------------
cat ("Data begin at line [default = 1]: ");
first_line_data = as.character(readLines("stdin", n=1));

if(first_line_data == '' | is.na(as.integer(first_line_data))){
    first_line_data = 1;
}
first_line_data = first_line_data - 1;

csv_data = "";
csv_data = csv_adjust(read.csv(file = path, sep = ";", skip = first_line_data)); #mudar o nome

cat ("Ploid Number: [default = 2]: ");
ploid = as.integer(readLines("stdin", n=1));

if(is.na(ploid)) {
    ploid = 2;
    print(paste0("Using ", ploid, "-ploid"));
}

loci_names = names(csv_data);
loci_names <- sort(loci_names[seq(1,length(loci_names), ploid)]); #Selects just one value for each locus(sorted)
total_loci = length(names(csv_data))/ploid;
total_individuals = nrow(csv_data);

print(paste("Your input data has", total_loci, "loci", "and", total_individuals, "individuals"));

cat ("Choose multiple number of sample size (n) for resamples [default = 5]: ");
n_minimum = as.character(readLines("stdin", n=1));

if(n_minimum == '') {
    n_minimum = 5;
    print("Using N = 5.");
}else{
    if(n_minimum > total_individuals) {
        print("Your N is higher than the total number of individuals in the original population!");
        print("Try Again.");
        quit();
    }else{
        if(n_minimum < 5){
            n_minimum = 5;
            print("N should not be less than 5!");
            print("Using N = 5.");
        }
    }
}

cat ("Number of resampling for each N class [default = 50]: ");
repeat_N = as.character(readLines("stdin", n=1));

if(repeat_N == '') {
    repeat_N = 50;
    print("Using 50 repeats.");
}else{
    if(repeat_N > 10000){
        print("Number of repetitions can not exceed 10.000!!!");
        quit();
    }
}

#frequence below which alleles are regarded as rare 
cat ("Minimal frequence of alleles to be preserved [default = 0.05]: "); 
freqMin = as.character(readLines("stdin", n=1));

if(freqMin == ''){
    freqMin = 0.05;
}else{
    if(freqMin > 0.90){
        freqMin = 0.90;
    }else{
        if(freqMin < 0.0001){
            freqMin = 0.0001;
        }
    }
}

print(paste("Using", freqMin, "as minimal frequence."));

if(os == "Linux") {
    system("clear");
}else{
    system("cls");
}

print(paste("Number for loci:", total_loci));
print(paste("Individuals:", total_individuals));

allele_rich <- gnot_numbers <- allele_numbers <- array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
names(allele_rich) <- names(gnot_numbers) <- names(allele_numbers) <- loci_names; #cria o nome da coluna

he_for_range = ho_for_range = array(rep(list(5), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
names(he_for_range) <- names(ho_for_range) <- loci_names;
# for(locus in loci_names) {
#         #he_for_range[[locus]] = ho_for_range[[locus]] = vector();
#         ho_for_range[[locus]] <- append(ho_for_range[[locus]], 5);
#         he_for_range[[locus]] <- append(he_for_range[[locus]], 5);
# }
start = Sys.time();
make_hash_count(csv_data); 
alleles_n = allele_count(allele_numbers);
end = Sys.time();
# print("Make_hash");
# print(end-start);
# quit();

print(paste("Number of alleles:", alleles_n));
print(paste("N class multiple:", n_minimum));
print(paste("Resample per N:", repeat_N));
print(paste("Alleles rares: bellow", freqMin));

cat ("Press any key to continue"); 
#readLines("stdin", n=1);
#print("--- Hashes ---");

allelic_richness(allele_numbers); #MUDAR PARA RETORNAR

orig_allele_numbers = allele_numbers;
orig_gnot_numbers = gnot_numbers;
orig_allele_rich = allele_rich;
#start = Sys.time();
orig_allele_freq = freq_calc(orig_allele_numbers);
orig_gnot_freq = freq_calc(orig_gnot_numbers);  #0.02588487 secs

#end = Sys.time();
#print("Tempo: ");
#print(end-start);
#print(orig_allele_freq);
#print(orig_gnot_freq);
#quit();
aux_array = freq_5(orig_allele_freq);
orig_five_percent_ref = aux_array[[1]]; five_percent_n = aux_array[[2]];

orig_five_percent = orig_five_percent_ref;

samples_five_percent = array(rep(list(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
names(samples_five_percent) <- loci_names;

freq_arrays_complete = list(); #VERIFICAR COMO ISSO VAI SER MANIPULADO DEPOIS

#print("Original frequencies:");
#print(orig_allele_freq);
#print(orig_gnot_freq);

aux_array = find_less_more(orig_allele_freq);
less_frq_all = aux_array[[1]]; more_frq_all = aux_array[[2]];
#print("Less and more:");
#print(paste(less_frq_all, "       ", more_frq_all));

he = ho = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
names(he) <- names(ho) <- loci_names;

# Takes the expct gnot freq for original data
expct_gnot_freq = expct_GNOT_freq(orig_allele_freq);
orig_expct_gnot_freq = expct_gnot_freq;

#print("Genotypes expected frequencies:");
#print(orig_expct_gnot_freq);

he_calc(orig_allele_freq, total_individuals, '1'); #MUDAR PARA RETORNAR
he_orig = he;
multi_he_orig = local_multi_he;

#FOR TEST ONLY
# print("He_ORIGINAL:");
# sorted_he_orig = sort(names(he_orig));
# for(locus in sorted_he_orig){
#     print(paste(locus, "=", he_orig[[locus]]));
# }
# print(paste("Multi_He_ORIGINAL =", multi_he_orig));

a = lapply(loci_names, function(locus){
    he[[locus]] = 0;
});

#-----------DO NOT SPLIT THIS CODE !!!-----------------
ho_calc(orig_gnot_freq); #MUDAR PARA RETORNAR
ho_range(ho); #MUDAR PARA RETORNAR
ho_orig = ho;
multi_ho_orig = local_multi_ho;
#------------------------------------------------------
#FOR TEST ONLY
# print("Ho_Original: ");
# sorted_ho_orig = sort(names(ho_orig));
# for(locus in sorted_ho_orig) {
#     print(paste(locus, "=", ho_orig[[locus]]));
# }
# print(paste("Multi_Ho_ORIGINAL =", multi_ho_orig));

print(" --- Initializing resampling !! ---");

stop_sign = 0;
repeats = 0;
allele_freq = array();
gnot_freq = array();
n_size = n_minimum;
#ht = hash();
multi_ht = 0;
multi_he = multi_ho = 0;
he_diff_final = ho_diff_final = vector();

he_diff_final <- append(he_diff_final, paste0(total_individuals, "\n", total_loci, "\n", n_minimum));
ho_diff_final <- append(ho_diff_final, paste0(total_individuals, "\n", total_loci, "\n", n_minimum));

he_mean = vector();
he_mean_final = vector(); # Mean, sd, min and max multilocus He for each n class
he_mean_final <- append(he_mean_final, paste0('#n', "\t", 'Mean_He', "\t", 'He_sd', "\t", 'min_He', "\t", 'max_He'));

ho_mean = vector();
ho_mean_final = vector(); # Mean, sd, min and max multilocus He for each n class
ho_mean_final <- append(ho_mean_final, paste0('#n', "\t", 'Mean_Ho', "\t", 'Ho_sd', "\t", 'min_Ho', "\t", 'max_Ho'));

mean_allele_freq = orig_allele_freq; #Copy orig_allele_freq to get structure of the matrix

global_freq_diffs = vector();
global_freq_diffs <- append(global_freq_diffs, paste0(total_individuals,"\n",alleles_n,"\n",n_minimum));

global_less_more_data = vector();
global_less_more_data <- append(global_less_more_data, paste0(total_individuals,"\n2\n",n_minimum));

global_fiverpercent_rate = vector();
global_fiverpercent_rate <- append(global_fiverpercent_rate, paste0(total_individuals,"\n",total_loci+1,"\n",n_minimum));

nei_fst = weir_theta = 0;
fst_array = vector(); # Fst values for resamples in a same N classes
fst_table = vector(); # Fst & SD for all N classes
neiD_array = vector(); # Nei distances for each N class
neiD_table = vector();
rogersD_array = vector();
rogersD_table = vector();
he_data = ho_data = vector(); # Mean, sd, min and max He/Ho for each locus in each n class
he_data <- append(he_data, paste0(total_individuals,"\n",total_loci,"\n",n_minimum));
ho_data <- append(ho_data, paste0(total_individuals,"\n",total_loci,"\n",n_minimum));

while(n_size <= total_individuals) {
    print(paste0("Calculating data for ", n_size, " Class."));
    start = Sys.time();
    he_mean = ho_mean = fst_array = neiD_array = rogersD_array = vector();
    repeats = repeat_N;
    
    five_percent_sample_count = ""; # string with number of presence of 5% freq alleles for each resample
    nClass_fiverpercent_count = vector(); #array to keep number of alleles w/ original freq >= 0.05
                                          #that are also present in each resemple for this N class
    
    #mean_allele_freq = hash();
    sorted_orig_allele_freq = sort(names(orig_allele_freq));

    for(locus in sorted_orig_allele_freq) {
        sorted_locus = sort(names(orig_allele_freq[[locus]]));

        for(allele in sorted_locus) {
            mean_allele_freq[[locus]][allele] = list(5);
        }
    }

    # for(locus in loci_names) {
    #     he_for_range[[locus]] = ho_for_range[[locus]] = array(5);
    # }
    he_for_range = ho_for_range = array(rep(list(5), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(he_for_range) <- names(ho_for_range) <- loci_names;

    #repeats = 2;
    
    freq_arrays_complete[[as.character(n_size)]] = list();
    start1 = Sys.time();
    while(repeats > 0) { # Do iterations for \$repeats resamples.
        allele_freq = gnot_freq = ht = vector();
        multi_ht = multi_he = multi_ho = nei_fst = weir_theta = 0;

        # for(locus in loci_names){
        #     he_diff[[locus]] = ho_diff[[locus]] = 5; # Defining structure of HoA 
        # }
        he_diff <- ho_diff <- array(rep(list(5), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
        names(he_diff) <- names(ho_diff) <- loci_names;

        # Resampling data !!
        sample1 = random_sampling(n_size);
        
        # Calculations for mixed data !!
        #mixed_data = vector(); # Merges original genotype table...
        mixed_data = rbind(sample1, csv_data); #VERIFICAR OS DADOS DISSO
        # if(n_size == 25){
        #     print(sample1);
        #     quit();
        #     print("que caraio");
        # }
        
        ht = mixed_calcs(mixed_data); # Calculates Ht;

        # Calculations for resample !!
        sample_calcs(sample1);
        # if(n_size == 30){
        #     #quit();
        # }
        #    print "\n \$Nei_Fst = $Nei_Fst \n"; # FOR TEST ONLY
        fst_array <- append(fst_array, nei_fst);

        #VERIFICAR ESSAS FUNÇÕES SE PRECISAM MESMO
        #sample1 = sample_for_R(sample1);
        # print(n_minimum);
        #freq_arrays_complete = array(rep(list(), 100/n_minimum-1, c(100/n_minimum)-1, 1)); #cria um array de tamanho loci
        #names(freq_arrays_complete) <- c(seq(n_size, 100-n_minimum, n_minimum));
        # print(freq_arrays_complete);
        # quit();
        
        r_input_HoH_HoA(allele_freq); #VERIFICAR ESSA FUNÇÃO PARA RETORNO

        samples_five_ref = 0;
        aux_array = count_freq_5_sample(allele_freq, samples_five_percent);
        
        five_percent_sample_count = aux_array[[1]]; samples_five_ref = aux_array[[2]];

        nClass_fiverpercent_count <- append(nClass_fiverpercent_count, five_percent_sample_count);
        samples_five_percent = samples_five_ref;
        
        repeats = repeats - 1;
    }
    end1 = Sys.time();
    #print(samples_five_percent);
    #aux_array = mean_and_sd(nClass_fiverpercent_count); #VERIFICAR ISSO, PROVAVELMETNE NÃO PRECISA DESSA FUNÇÃO
    nClass_fiverpercent_mean = sprintf("%.5f", mean(nClass_fiverpercent_count));
    nClass_fiverpercent_sd = sprintf("%.5f", sd(nClass_fiverpercent_count));

    global_fiverpercent_rate <- append(global_fiverpercent_rate, paste0(n_size, "\t", "ML", "\t", nClass_fiverpercent_mean, "\t", nClass_fiverpercent_sd));
    #print(global_fiverpercent_rate);
    
    
    sorted_samples = sort(names(samples_five_percent));
    for(locus in sorted_samples) {

        nClass_fiverpercent_mean = sprintf("%.5f", mean(samples_five_percent[[locus]]));
        nClass_fiverpercent_sd = sprintf("%.5f", sd(samples_five_percent[[locus]]));
        
        global_fiverpercent_rate <- append(global_fiverpercent_rate, paste0(n_size, "\t", locus, "\t", nClass_fiverpercent_mean, "\t", nClass_fiverpercent_sd));
        samples_five_percent[[locus]] = vector(); # Clean \%SAMPLES_five_percent values for use in next N class
    }

    local_freq_diffs = freq_diff_calc(mean_allele_freq);
    #print(local_freq_diffs);
    
    global_freq_diffs <- append(global_freq_diffs, local_freq_diffs);

    local_less_more_data = less_more_data(mean_allele_freq);
    global_less_more_data <- append(global_less_more_data, local_less_more_data);

    he_mean = sort(he_mean);
    mean_he = sprintf("%.5f", mean(he_mean)); sd_he = sprintf("%.5f", sd(he_mean));
    
    he_mean_final <- append(he_mean_final, paste0(n_size, "\t", mean_he, "\t", sd_he, "\t", he_mean[1], "\t", he_mean[length(he_mean)]));
    he_mean = vector();

    sorted_he_range = sort(names(he_for_range));
    for(locus in sorted_he_range){
        #he_for_range[[locus]] = he_for_range[[locus]][-1];
        local_array = sort(he_for_range[[locus]][-1]);
        #he_med = mean(he_for_range[[locus]]); he_sd = sd(he_for_range[[locus]]);
        he_med = sprintf("%.5f", mean(local_array)); he_sd = sprintf("%.5f", sd(local_array));

        #local_array = sort(he_for_range[[locus]]); #VERIFICAR SE É SOMENTE NAMES OU TUDO
        he_min = sprintf("%.5f", local_array[1]);
        he_max = sprintf("%.5f", local_array[length(local_array)]);

        he_data <- append(he_data, paste0(n_size, "\t", locus, "\t", he_med, "\t", he_sd, "\t", he_min, "\t", he_max));
        
        #he_data <- append(he_data, paste(he_min, "\t", he_max, "\n"));
    }
    
    ho_mean = sort(ho_mean);
    mean_ho = sprintf("%.5f", mean(ho_mean)); sd_ho = sprintf("%.5f", sd(ho_mean));
    ho_mean_final <- append(ho_mean_final, paste0(n_size, "\t", mean_ho, "\t", sd_ho, "\t", ho_mean[1], "\t", ho_mean[length(ho_mean)]));
    ho_mean = vector();
    #ho_mean_final <- append(ho_mean_final, c(ho_mean[1], "\t")); #Extract min he
    #ho_mean_final <- append(ho_mean_final, c(ho_mean[length(ho_mean)], "\n")); #Extract max He VERIFICAR ESSES DOIS


    sorted_ho_range = sort(names(ho_for_range));
    for(locus in sorted_ho_range){
        #he_for_range[[locus]] = he_for_range[[locus]][-1];
        local_array = sort(ho_for_range[[locus]][-1]);
        #he_med = mean(he_for_range[[locus]]); he_sd = sd(he_for_range[[locus]]);
        ho_med = sprintf("%.5f", mean(local_array)); ho_sd = sprintf("%.5f", sd(local_array));

        #local_array = sort(he_for_range[[locus]]); #VERIFICAR SE É SOMENTE NAMES OU TUDO
        ho_min = sprintf("%.5f", local_array[1]);
        ho_max = sprintf("%.5f", local_array[length(local_array)]);

        ho_data <- append(ho_data, paste0(n_size, "\t", locus, "\t", ho_med, "\t", ho_sd, "\t", ho_min, "\t", ho_max));
        
        #he_data <- append(he_data, paste(he_min, "\t", he_max, "\n"));
    }
    #print(ho_data);

    he_diff_final <- append(he_diff_final, het_diff_push(he_diff));
    ho_diff_final <- append(ho_diff_final, het_diff_push(ho_diff));

    mean_fst = sprintf("%.5f", mean(fst_array)); sd_fst = sprintf("%.5f", sd(fst_array));
    fst_table <- append(fst_table, paste0(n_size, "\t", mean_fst, "\t", sd_fst));
    #print(paste("Fst_table:", fst_table));

    mean_neiD = sprintf("%.5f", mean(neiD_array)); sd_neiD = sd(neiD_array);
    neiD_table <- append(neiD_table, paste0(n_size, "\t", mean_neiD, "\t", ifelse(is.na(sd_neiD), 0, sprintf("%.5f", sd_neiD))));
    #print(paste("neiD_table:", neiD_table));

    mean_rogersD = sprintf("%.5f", mean(rogersD_array)); sd_rogersD = sd(rogersD_array);
    rogersD_table <- append(rogersD_table, paste0(n_size, "\t", mean_rogersD, "\t", ifelse(is.na(sd_rogersD), 0, sprintf("%.5f", sd_rogersD))));
    #print(paste("RogersD_table:", rogersD_table));

    #print("He_mean_final");
    #print(he_mean_final);

    #print("Ho_mean_final");
    #print(ho_mean_final);
    print("Printing results into output files...");
    #----SAVING some RESULTS TO OUTPUT FILES-----#
    save_array_output(global_fiverpercent_rate, '1-5percent_rate');
    save_array_output(global_freq_diffs, '2-freq_dif');
    save_array_output(global_less_more_data, '3-freq_impact');
    save_array_output(he_data, '4-He_impact');
    save_array_output(ho_data, '10-Ho_impact');
    save_array_output(he_diff_final, '12-He_diff');
    save_array_output(he_mean_final, '5-meanHe_impact');
    save_array_output(ho_mean_final, '9-meanHo_impact');
    save_array_output(ho_diff_final, '11-Ho_diff');
    save_array_output(fst_table, '6-Fst');
    save_array_output(neiD_table, '7-Nei');
    save_array_output(rogersD_table, '8-Roger');

    #-----RESULTS SAVED---------------------#

    #freq_arrays_forR = r_input_HoHoA_array(freq_arrays_complete);

    #----SAVING r RESULTS TO OUTPUT FILES-----#
    #save_array_output(freq_arrays_forR, 'ALL_freq_R_fig3');
    #-----RESULTS SAVED---------------------#

    #---------------------------------------------------------------
    n_size = n_size + n_minimum; # Be carefull if move this code!!!

    if (stop_sign != 0){
        #make_temp_files();
        plot_graphs();
        print("Ploting graphs...");
        print("Finish !");
        quit();
    }else {
        if (n_size >= total_individuals){
            n_size = total_individuals;
            stop_sign = 2;
        }
    }
    #----------------------------------------------------------------
    end = Sys.time();
    tempoTotal = end-start;
    tempoInterno = end1-start1;
    print(paste0("Tempo total: ", tempoTotal));
    print(paste0("Tempo while interno: ", tempoInterno));
    print(paste0("Tempo sem while interno: ", tempoTotal-tempoInterno));
    quit();
}
