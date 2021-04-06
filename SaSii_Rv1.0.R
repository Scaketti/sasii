#-------REQUESTS-------------------------
if (!require(Rcpp)) install.packages('Rcpp')
if (!require(readr)) install.packages('readr')
#if (!require(microbenchmark)) install.packages('microbenchmark')
#----------------------------------------
# ------SOME GLOBAL VARIABLES------------
args <- commandArgs(trailingOnly = TRUE);

os = Sys.info()['sysname'];
path = '';
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
# Creates output file with results from script
save_array_output <- function(lines, mid_name){
    output_name = paste0(mid_name,'.sso');
    path = "results/";
    file.create(paste0(path, output_name));

    fileConn<-file(paste0(path, output_name));
    writeLines(lines, fileConn);
    close(fileConn);
}

#----------------------------------------
# Make adjustment in csv data
# For type 1 csv file, just changes second locus column name to {{locus_name}}.2 to aid in code optmization
# For type 2 csv file, create a new dataframe to change the individual allele order, just for padronization
csv_adjust <- function(csv_data){
    locus = "";
    #csv_data = csv_data[-length(csv_data)];

    #Remove undesirable columns
    if(first_column_data > 1){
        csv_data = csv_data[,-1:-(first_column_data-1)];
    }

    #If csv individual alleles is in just one line, copy the name of every locus adding ".2"
    if(str_type == 1){
        for(lines in names(csv_data)){
            if(grepl('X', lines, fixed = TRUE)) {
                names(csv_data)[names(csv_data) == lines] <- paste0(locus, '.2');
            }else{locus = lines;}
        }
    }else{
        csv_data[csv_data==indef] = "indef"; #Look for each indefined allele and changes its value
        new_data = data.frame(matrix(ncol = ncol(csv_data) * 2, nrow = 0)); #Creates a temp variable to change csv style
        dplc_locus = array(unlist(rbind(names(csv_data),paste0(names(csv_data), ".2")))); #Duplicate locus columns
        names(new_data) <- dplc_locus;

        #For each individual, change alleles order into one line pattern
        lapply(seq(1,nrow(csv_data), 2), function(indv){
            new_indv = data.frame(matrix(data = unlist(rbind(csv_data[indv,],csv_data[indv+1,])), nrow = 1, ncol = ncol(csv_data)*2)); #Create a "list" of alternating values between line A and B
            names(new_indv) <- dplc_locus;
            new_data <<- rbind(new_data, new_indv); #Bind new generated row into the new csv_data
        });

        csv_data = new_data;
    }
    return (csv_data);
}
#----------------------------------------
#Create two matrixes with alleles and genotypes counted values
count_allele_gnot <- function(data){
    a = lapply(loci_names, function(locus){
        alleles = unlist(subset(data, TRUE, c(locus, paste0(locus,".2"))));

        #Verify csv file to remove indefined alleles from counting
        if(str_type == 2){
            alleles = alleles[alleles != "indef"];
            counted_alleles = table(droplevels(alleles[alleles != "indef"]));
        }else{
            counted_alleles = table(alleles);
        }

        allele_numbers[[locus]] <<- list(counted_alleles[is.numeric(counted_alleles)]); 

        counted_gnot = table(paste0(alleles[seq(1, length(alleles)/2)], "-", alleles[seq(length(alleles)/2 + 1, length(alleles))]))
        gnot_numbers[[locus]] <<- list(counted_gnot);
    });
}
#----------------------------------------
# Count alleles in population
allele_count <- function(allele_hash){
    count = 0;
    a = lapply(allele_hash, function(locus){
        lapply(names(locus[[1]]), function(allele){
            if(allele != '0'){
                count <<- count + 1;
            }
        })
    });
    return(count);
}
#----------------------------------------
#Counts alleles numbers for each locus to obtain allelic richness
allelic_richness <- function(allele_hash){
    a = lapply(names(allele_hash), function(locus){
        allele_rich[[locus]] <<- length(allele_hash[[locus]]);
    });
}
#----------------------------------------
# Old freq_calc function
# freq_calc <- function(numbers_hash){
#     freq_hash = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
#     names(freq_hash) <- loci_names;
#     sorted_numbers_hash = sort(names(numbers_hash));

#     a = lapply(sorted_numbers_hash, function(locus){
#         sorted_locus_hash = sort(names(numbers_hash[[locus]][[1]]));
#         #Alleles is a list without lost data
#         good_data = as.character(unlist(
#             lapply(sorted_locus_hash, 
#                    function(x) x[(which(x != '000' & x != '000-000' & 
#                                         x != '0' & x != '0-0'))])));

#         freq_hash[[locus]] <<- list(numbers_hash[[locus]][[1]][good_data]); #Creates a copy of the variable type
#         total_alleles = sum(numbers_hash[[locus]][[1]][good_data]);

#         freq_hash[[locus]] <<- numbers_hash[[locus]][[1]][good_data] / total_alleles; #Calculate the frequency for each element
#     });
#     return(freq_hash);
# }
#C_FREQ_CALC
# Calculate frequencies for alleles and genotypes. Get counted alleles and genotypes divide by the sum of the counted alleles in locus
# Ex: One locus has two alleles: allele_158 -> 2 and allele_202 -> 8. The sum of alleles is 10. 
# So the respectives frequencies are: allele_158 -> 0.2 and allele_202 -> 0.8.
cppFunction('
    List freq_calc(List numbers_hash){
        CharacterVector locus_name(numbers_hash.names());
        List freq_hash = List::create();
        int total_alleles = 0;
        for(auto pLocus : locus_name){
            std::string locus = as<std::string>(pLocus);
            List alleles = as<List> (numbers_hash[locus])[0];
            CharacterVector alleles_name(alleles.names());
            NumericVector good_alleles = NumericVector::create();

            for(auto pAllele : alleles_name){
                std::string allele = as<std::string>(pAllele);
                if(allele != "0" && allele != "0-0") good_alleles.push_back(alleles[allele], allele);
            }
            total_alleles = sum(good_alleles);
            NumericVector res = good_alleles / total_alleles;
            res.names() = good_alleles.names();

            freq_hash.push_back(res, locus);
            
        }
        return freq_hash;
    }
')
#----------------------------------------
# Get smaller and bigger frequencies values from population, extrating their respectives locus and allele.
find_less_more <- function(freqs){
    #Extract the position of the less freq allele
    less_pos = match(min(unlist(freqs)), unlist(freqs))
    less_locus = strsplit(names(unlist(freqs)[less_pos]), "\\.")[[1]][1];
    less_allele = strsplit(names(unlist(freqs)[less_pos]), "\\.")[[1]][2];

    #Extract the position of the more freq allele
    more_pos = match(max(unlist(freqs)), unlist(freqs))
    more_locus = strsplit(names(unlist(freqs)[more_pos]), "\\.")[[1]][1];
    more_allele = strsplit(names(unlist(freqs)[more_pos]), "\\.")[[1]][2];

    local_less = paste0(less_locus, "-", less_allele);
    local_more = paste0(more_locus, "-", more_allele);
    return (list(local_less, local_more));
}
#----------------------------------------
# Sort frequencies to look up smaller and bigger frequencies for resampling step
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
# Calculate expected genotypes frquencies from HoH with allelic freq
expct_GNOT_freq <- function(alelic_freq){ 
    EXPCT_GNOT_freq = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1));
    names(EXPCT_GNOT_freq) <- loci_names; # Hash with expected genotypes frequencies

    sorted_alelic_freq = sort(names(alelic_freq));

    a = lapply(sorted_alelic_freq, function(locus){
        exp_gnot_freq = 0;
        generated_genotypes = generate_genotypes(names(alelic_freq[[locus]]));

        EXPCT_GNOT_freq[locus] <<- list(array(rep(array(), length(generated_genotypes)), c(length(generated_genotypes), 1)));
        names(EXPCT_GNOT_freq[[locus]]) <<- generated_genotypes;

        lapply(generated_genotypes, function(genotype){
            
            aux_array = strsplit(genotype, '-')[[1]];
            allele1 = aux_array[1]; allele2 = aux_array[2];

            if(allele1 == allele2){
                exp_gnot_freq = alelic_freq[[locus]][allele1] * alelic_freq[[locus]][allele2];
            }else {
                if (allele1 != allele2) {
                    exp_gnot_freq = 2 * (alelic_freq[[locus]][allele1] * alelic_freq[[locus]][allele2]);
                }
            }

            EXPCT_GNOT_freq[[locus]][genotype] <<- exp_gnot_freq;
        });
    });
    return(EXPCT_GNOT_freq);
}
#----------------------------------------
# Create list of possible genotypes from an hash: allele => value
generate_genotypes <- function(allelic_freq_inlocus) {
    generated_genotypes = array(); # Array to save genotypes list

    sorted_allelic_freq_inlocus = sort(allelic_freq_inlocus);
    a = lapply(sorted_allelic_freq_inlocus, function(allele1){
        lapply(sorted_allelic_freq_inlocus, function(allele2) {
            genotype = paste0(allele1, '-', allele2);

            if(!is.na(generated_genotypes[1])){
                generated_genotypes <<- append(generated_genotypes, genotype);
            }else{
                generated_genotypes <<- array(genotype);
            }
        });
        sorted_allelic_freq_inlocus <<- sorted_allelic_freq_inlocus[-1];
    });
    return(generated_genotypes);
}
#----------------------------------------
# Find observed heterosygosity for data frame
ho_calc <- function(gnot_freq) { 
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
    });
    sorted_ho = sort(names(ho));

    ho_array_local = vector();
    a = lapply(sorted_ho, function(locus){
        ho_array_local <<- append(ho_array_local, ho[[locus]]);
    });
    
    mean = sprintf("%.5f", mean(ho_array_local)); SD = sprintf("%.5f", sd(ho_array_local));

    local_multi_ho <<- mean;
}
#----------------------------------------
# Get random individuals from population for sampling
random_sampling <- function(n) {
    count = n; 
    total_lines = total_individuals;
    sample = vector();
    
    sorted_numbers = sample(1:nrow(csv_data), n);
    #sorted_numbers = rev(1:n)#sample(rev(1:n), n);
    sample = csv_data[sorted_numbers,];

    return(sample);
}
#----------------------------------------
# Calculate important information for each sample
sample_calcs <- function(sample) {
    allele_rich <<- gnot_numbers <<- allele_numbers <<- array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(allele_rich) <<- names(gnot_numbers) <<- names(allele_numbers) <<- loci_names; #cria o nome da coluna

    count_allele_gnot(sample);
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

    ho_calc(gnot_freq);
    ho_range(ho);
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

    count_allele_gnot(sample);
    allele_freq <<- freq_calc(allele_numbers);
    gnot_freq = freq_calc(gnot_numbers);
    EXPCT_GNOT_freq = expct_GNOT_freq(allele_freq);
    he_calc(allele_freq, n_size, '1');
    mix_he = he;
    multi_ht <<- local_multi_he;

    return(mix_he);
}
#----------------------------------------
#Needs a HoH with allelic frequencies
he_calc <- function(local_hash, number_id, type) { 
    somatory = 0;
    r = 0; #Number o loci

    sorted_local_hash = sort(names(local_hash));

    a = lapply(sorted_local_hash, function(locus){
        somatory = sum(local_hash[[locus]] ** 2);

        he_simple = 1 - somatory;
        he[[locus]] <<- he_simple * ((2*number_id)/((2*number_id)-1));
    
        if(type == '2') {
            he_for_range[[locus]] <<- append(he_for_range[[locus]], he[[locus]]); 
        }
    });

    local_multi_he <<- 0;
    he_array_local = vector();
    sorted_he = sort(names(he));

    a = lapply(sorted_he, function(locus){
        he_array_local <<- append(he_array_local, he[[locus]]);
    });

    mean = sprintf("%.5f", mean(he_array_local)); SD = sprintf("%.5f", sd(he_array_local));
    local_multi_he <<- mean;
}
#----------------------------------------
#Calcs Nei Fst (1973). Needs $Multi_Ht, $Multi_He_ORIGINAL, $Multi_He, $n_size.
calc_nei_fst <- function(ht, hsA, hsB, nB) { 
    nA = total_individuals;
    nei_fst = 0;
    
    upper1a = ((as.numeric(nA) * as.numeric(hsA)) + (as.numeric(nB) * as.numeric(hsB))) / (as.numeric(nA) + as.numeric(nB));
    upper1b = as.numeric(ht) - upper1a;
    nei_fst = upper1b / as.numeric(ht);

    return (nei_fst);
}
#----------------------------------------
#Nei distance (1972)
nei_d <- function(freq2) { 
    freq1 = orig_allele_freq;

    jx = jy = jxy = l = 0; 
    jx_array = jy_array = jxy_array = vector();

    sorted_freq1 = sort(names(freq1));
    a = lapply(sorted_freq1, function(locus) {
        l <<- l + 1;

        sorted_locus = sort(names(freq1[[locus]]));
        lapply(sorted_locus, function(allele) {
            if ( is.na(freq2[[locus]][allele])){
                freq2[[locus]][allele] <<- 0;
            }
            jx <<- freq1[[locus]][allele] ** 2;
            jx_array <<- append(jx_array, jx);
            jy <<- freq2[[locus]][allele] ** 2;
            jy_array <<- append(jy_array, jy);
        });
    });

    a = lapply(sorted_freq1, function(locus) {
        sorted_locus = sort(names(freq1[[locus]]));

        lapply(sorted_locus, function(allele) {
            if(is.na(freq2[[locus]][allele]) || freq2[[locus]][allele] == 0){ jxy <<- 0;}
            else {
                jxy <<- freq1[[locus]][allele] * freq2[[locus]][allele];
            }
            jxy_array <<- append(jxy_array, jxy);
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
# Calculares modified Rogers distance. Needs HoH of resample allele frequencies. Roger (1972), Wright (1978), GOodman and Stuber (1983) 
rogers_d <- function(freq2) { 

    freq1 = orig_allele_freq;

    diff_sum = 0; l = 0;
    diff_sum_array = vector();

    sorted_freq1 = sort(names(freq1));

    a = lapply(sorted_freq1, function(locus){
        l <<- l + 1;
        diff_array = vector();
        #print(locus);
        sorted_locus = sort(names(freq1[[locus]]));

        lapply(sorted_locus, function(allele) {
            #if(is.na(freq2[[locus]][allele]) || freq2[[locus]][allele] == 0) freq2[[locus]][allele] <<- 0; #Treat cases where allele do not exist
            diff = freq1[[locus]][allele] - freq2[[locus]][allele];
            diff = diff ** 2;
            
            diff_array <<- append(diff_array, diff);
        });

        diff_sum = sum(diff_array);
        diff_sum_array <<- append(diff_sum_array, diff_sum);
    });
    
    rogersD = sqrt(sum(diff_sum_array));
    rogersD = rogersD / sqrt(2 * l);

    return (rogersD);
}
#----------------------------------------
mean_freq_prep1 <- function(freqs) {
    sorted_orig_allele_freq = sort(names(orig_allele_freq));

    a = lapply(sorted_orig_allele_freq, function(locus){
        sorted_locus = sort(names(orig_allele_freq[[locus]]));

        lapply(sorted_locus, function(allele){
            if(is.na(allele_freq[[locus]][allele])){
                allele_freq[[locus]][allele] <<- 0;
            }
            mean_allele_freq[[locus]][allele] <<- list(append(unlist(mean_allele_freq[[locus]][allele]), allele_freq[[locus]][allele]));
        });
    });
}
#----------------------------------------
freq_diff_calc <- function(freqs) {
    mean_freqs = freqs;
    sorted_freqs = sort(names(freqs));
    freqs_differences = vector();

    for(locus in sorted_freqs) {
        sorted_locus = sort(names(freqs[[locus]]));

        for(allele in sorted_locus) {
            mean_freqs[[locus]][allele] = sprintf('%.5f', sum(freqs[[locus]][[allele]][-1])/length(freqs[[locus]][[allele]][-1]));

            diff = orig_allele_freq[[locus]][allele] - as.double(mean_freqs[[locus]][[allele]]);
            diff = abs(diff);

            if(allele < 100) {allele = paste0('_', allele);}

            freqs_differences <- append(freqs_differences, paste0(n_size, "\t", paste0(locus, ".", allele), "\t", sprintf('%.5f', diff)));
        }
    }
    return (freqs_differences);
}
#----------------------------------------
freq_5 <- function(freq_hash) {
    five_percent_hash = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1));
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
    
    a = lapply(sorted_freq_hash, function(locus) {
        local_array = as.array(unlist(reference_hash[[locus]]));
        sorted_freq_locus = sort(names(freq_hash[[locus]]));
        local_count = 0;

        lapply(sorted_freq_locus, function(key){
            #print(paste(key, "= ", freq_hash[[locus]][key], key %in% local_array));
            if((freq_hash[[locus]][key] != 0) && (key %in% local_array)) { #verificar se está correto
                local_count <<- local_count + 1;
                general_local_count <<- general_local_count + 1;
            }
        });
        local_count = local_count / length(local_array);
        sample_local[[locus]] <<- append(sample_local[[locus]], local_count);
    })
    general_local_count = general_local_count / local_ref_count;
    
    return (list(general_local_count, sample_local));
}
#----------------------------------------
r_input_HoHoA_array <- function(target_hash1) {
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
plot_graphs <- function() {
    source("SaSii_plots.R");
}
#----------------------------------------
het_diff_calc <- function(hash, local_hash, orig_hash) {
    sorted_orig_hash = sort(names(orig_hash));
    for(locus in sorted_orig_hash) {
        local_hash[[locus]] <- append(local_hash[[locus]], abs(orig_hash[[locus]] - hash[[locus]]));
    }
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

if(file.exists("config")){
    config = read_lines("config", n_max = 10);
}else{
    config = NA;
}

if(!is.na(config[1])){
    path = strsplit(config[1], " ")[[1]][2];
    str_type = as.integer(strsplit(config[2], " ")[[1]][2]);
    indef = strsplit(config[3], " ")[[1]][2];
    first_line_data = as.integer(strsplit(config[4], " ")[[1]][2]);
    n_minimum = as.integer(strsplit(config[5], " ")[[1]][2]);
    repeat_N = as.integer(strsplit(config[6], " ")[[1]][2]);
    freqMin = as.numeric(strsplit(config[7], " ")[[1]][2]);
    separator = as.character(strsplit(config[8], " ")[[1]][2]);
    nLines = as.numeric(strsplit(config[9], " ")[[1]][2]);
    first_column_data = as.integer(strsplit(config[10], " ")[[1]][2]);
    print(paste(separator, intToUtf8(separator)));
}

if(is.na(config[1])){
    if(length(args) == 0){
        cat("Put you input file path: ");
        path = as.character(readLines("stdin", n=1));
    }else{
        path = args;
    }
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
if(is.na(config[1])){
    cat ("Data begin at line [default = 1]: ");
    first_line_data = as.character(readLines("stdin", n=1));
}

if(first_line_data == '' | is.na(as.integer(first_line_data))){
    first_line_data = 1;
}
first_line_data = first_line_data - 1;

csv_data = "";

csv_data = csv_adjust(read.csv(file = path, sep = intToUtf8(separator), nrows = nLines, skip = first_line_data)); #mudar o nome

#if(is.na(ploid)) {
ploid = 2; #May be a future change for multiploid species
    #print(paste0("Using ", ploid, "-ploid"));
#}
loci_names = names(csv_data);
loci_names <- sort(loci_names[seq(1,length(loci_names), ploid)]); #Selects just one value for each locus(sorted)
total_loci = length(names(csv_data))/ploid;
total_individuals = nrow(csv_data);

print(paste("Your input data has", total_loci, "loci", "and", total_individuals, "individuals"));

if(is.na(config[1])){
    cat ("Choose multiple number of sample size (n) for resamples [default = 5]: ");
    n_minimum = as.character(readLines("stdin", n=1));
}

if(n_minimum == '') {
    n_minimum = 5;
    print("Using N = 5.");
}else{
    n_minimum = as.integer(n_minimum);
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

if(is.na(config[1])){
    cat ("Number of resampling for each N class [default = 50]: ");
    repeat_N = as.character(readLines("stdin", n=1));
}

if(repeat_N == '') {
    repeat_N = 50;
    print("Using 50 repeats.");
}else{
    repeat_N = as.integer(repeat_N); 
    if(repeat_N > 10000){
        print("Number of repetitions can not exceed 10.000!!!");
        quit();
    }
}

#frequence below which alleles are regarded as rare 
if(is.na(config[1])){
    cat ("Minimal frequence of alleles to be preserved [default = 0.05]: "); 
    freqMin = as.character(readLines("stdin", n=1));
}

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

count_allele_gnot(csv_data); 
alleles_n = allele_count(allele_numbers);

print(paste("Number of alleles:", alleles_n));
print(paste("N class multiple:", n_minimum));
print(paste("Resample per N:", repeat_N));
print(paste("Alleles rares: bellow", freqMin));

cat ("Press any key to continue\n"); 
readLines("stdin", n=1);
#print("--- Hashes ---");

allelic_richness(allele_numbers);

orig_allele_numbers = allele_numbers;
orig_gnot_numbers = gnot_numbers;
orig_allele_rich = allele_rich;

orig_allele_freq = freq_calc(orig_allele_numbers);
orig_gnot_freq = freq_calc(orig_gnot_numbers)

aux_array = freq_5(orig_allele_freq);
orig_five_percent_ref = aux_array[[1]]; five_percent_n = aux_array[[2]];

orig_five_percent = orig_five_percent_ref;

samples_five_percent = array(rep(list(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1));
names(samples_five_percent) <- loci_names;

freq_arrays_complete = list();

aux_array = find_less_more(orig_allele_freq);
less_frq_all = aux_array[[1]]; more_frq_all = aux_array[[2]];

he = ho = array(rep(array(), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1));
names(he) <- names(ho) <- loci_names;

# Takes the expct gnot freq for original data
expct_gnot_freq = expct_GNOT_freq(orig_allele_freq);
orig_expct_gnot_freq = expct_gnot_freq;

he_calc(orig_allele_freq, total_individuals, '1');
he_orig = he;
multi_he_orig = local_multi_he;

a = lapply(loci_names, function(locus){
    he[[locus]] <<- 0;
});

#-----------DO NOT SPLIT THIS CODE !!!-----------------
ho_calc(orig_gnot_freq);
ho_range(ho);
ho_orig = ho;
multi_ho_orig = local_multi_ho;
#------------------------------------------------------

print(" --- Initializing resampling !! ---");
stop_sign = repeats = multi_ht = multi_he = multi_ho = 0;
allele_freq = gnot_freq = array();
n_size = n_minimum;
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
    
    he_mean = ho_mean = fst_array = neiD_array = rogersD_array = vector();
    repeats = repeat_N;
    
    five_percent_sample_count = ""; # string with number of presence of 5% freq alleles for each resample
    nClass_fiverpercent_count = vector(); #array to keep number of alleles w/ original freq >= 0.05
                                          #that are also present in each resemple for this N class

    sorted_orig_allele_freq = sort(names(orig_allele_freq));

    a = lapply(names(sorted_orig_allele_freq), function(locus){
        sorted_locus = sort(names(orig_allele_freq[[locus]]));

        lapply(sorted_locus, function(allele){
            mean_allele_freq[[locus]][allele] <<- list(5);
        });
    });

    he_for_range = ho_for_range = array(rep(list(5), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
    names(he_for_range) <- names(ho_for_range) <- loci_names;
    freq_arrays_complete[[as.character(n_size)]] = list();

    while(repeats > 0) { # Do iterations for \$repeats resamples.
        # print(repeats);
        # print(" ")
        # print(" ")
        # print(" ")

        allele_freq = gnot_freq = ht = vector();
        multi_ht = multi_he = multi_ho = nei_fst = weir_theta = 0;

        he_diff <- ho_diff <- array(rep(list(5), ncol(csv_data)/ploid), c(ncol(csv_data)/ploid, 1)); #cria um array de tamanho loci
        names(he_diff) <- names(ho_diff) <- loci_names;

        # Resampling data !!
        sample1 = random_sampling(n_size);
        
        # Calculations for mixed data !!
        mixed_data = rbind(sample1, csv_data);
        ht = mixed_calcs(mixed_data); # Calculates Ht;
        sample_calcs(sample1);
        
        fst_array <- append(fst_array, nei_fst);
        samples_five_ref = 0;
        
        aux_array = count_freq_5_sample(allele_freq, samples_five_percent);
        
        five_percent_sample_count = aux_array[[1]]; samples_five_ref = aux_array[[2]];

        nClass_fiverpercent_count <- append(nClass_fiverpercent_count, five_percent_sample_count);
        samples_five_percent = samples_five_ref;
        
        repeats = repeats - 1;
    }
    #quit();
    nClass_fiverpercent_mean = sprintf("%.5f", mean(nClass_fiverpercent_count));
    nClass_fiverpercent_sd = sprintf("%.5f", sd(nClass_fiverpercent_count));

    global_fiverpercent_rate <- append(global_fiverpercent_rate, paste0(n_size, "\t", "ML", "\t", nClass_fiverpercent_mean, "\t", nClass_fiverpercent_sd));

    sorted_samples = sort(names(samples_five_percent));
    a = lapply(sorted_samples, function(locus){
        nClass_fiverpercent_mean = sprintf("%.5f", mean(samples_five_percent[[locus]]));
        nClass_fiverpercent_sd = sprintf("%.5f", sd(samples_five_percent[[locus]]));
        
        global_fiverpercent_rate <<- append(global_fiverpercent_rate, paste0(n_size, "\t", locus, "\t", nClass_fiverpercent_mean, "\t", nClass_fiverpercent_sd));
        samples_five_percent[[locus]] <<- vector(); # Clean \%SAMPLES_five_percent values for use in next N class
    });

    local_freq_diffs = freq_diff_calc(mean_allele_freq);
    
    global_freq_diffs <- append(global_freq_diffs, local_freq_diffs);

    local_less_more_data = less_more_data(mean_allele_freq);
    global_less_more_data <- append(global_less_more_data, local_less_more_data);

    he_mean = sort(he_mean);
    mean_he = sprintf("%.5f", mean(he_mean)); sd_he = sprintf("%.5f", sd(he_mean));
    
    he_mean_final <- append(he_mean_final, paste0(n_size, "\t", mean_he, "\t", sd_he, "\t", he_mean[1], "\t", he_mean[length(he_mean)]));
    he_mean = vector();

    sorted_he_range = sort(names(he_for_range));
    a = lapply(sorted_he_range, function(locus){
        local_array = sort(he_for_range[[locus]][-1]);
        he_med = sprintf("%.5f", mean(local_array)); he_sd = sprintf("%.5f", sd(local_array));
        he_min = sprintf("%.5f", local_array[1]);
        he_max = sprintf("%.5f", local_array[length(local_array)]);
        he_data <<- append(he_data, paste0(n_size, "\t", locus, "\t", he_med, "\t", he_sd, "\t", he_min, "\t", he_max));
    });
    
    ho_mean = sort(ho_mean);
    mean_ho = sprintf("%.5f", mean(ho_mean)); sd_ho = sprintf("%.5f", sd(ho_mean));
    ho_mean_final <- append(ho_mean_final, paste0(n_size, "\t", mean_ho, "\t", sd_ho, "\t", ho_mean[1], "\t", ho_mean[length(ho_mean)]));
    ho_mean = vector();

    sorted_ho_range = sort(names(ho_for_range));
    a = lapply(sorted_ho_range, function(locus){
        local_array = sort(ho_for_range[[locus]][-1]);
        ho_med = sprintf("%.5f", mean(local_array)); ho_sd = sprintf("%.5f", sd(local_array));
        ho_min = sprintf("%.5f", local_array[1]);
        ho_max = sprintf("%.5f", local_array[length(local_array)]);
        ho_data <<- append(ho_data, paste0(n_size, "\t", locus, "\t", ho_med, "\t", ho_sd, "\t", ho_min, "\t", ho_max));
    });

    he_diff_final <- append(he_diff_final, het_diff_push(he_diff));
    ho_diff_final <- append(ho_diff_final, het_diff_push(ho_diff));

    mean_fst = sprintf("%.5f", mean(fst_array)); sd_fst = sprintf("%.5f", sd(fst_array));
    fst_table <- append(fst_table, paste0(n_size, "\t", mean_fst, "\t", sd_fst));

    mean_neiD = sprintf("%.5f", mean(neiD_array)); sd_neiD = sd(neiD_array);
    neiD_table <- append(neiD_table, paste0(n_size, "\t", mean_neiD, "\t", ifelse(is.na(sd_neiD), 0, sprintf("%.5f", sd_neiD))));

    mean_rogersD = sprintf("%.5f", mean(rogersD_array)); sd_rogersD = sd(rogersD_array);
    rogersD_table <- append(rogersD_table, paste0(n_size, "\t", mean_rogersD, "\t", ifelse(is.na(sd_rogersD), 0, sprintf("%.5f", sd_rogersD))));

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
}
