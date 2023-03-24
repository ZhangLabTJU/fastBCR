df2fas = function(data){
  ids = data$ids
  seqs = data$germline
  res = paste0(">",ids,"\n",seqs)
  # write.table(res,
  #             file = file,
  #             row.names = F,
  #             col.names = F,
  #             quote = F)
  return(res)
}

#' Function: Generation of a fasta file consisting of germline sequences
#'
#' @description Randomly generate germline sequences and convert them to a fasta file.
#'
#' @param family_num The number of simulated clonal families.
#'
#' @return A faste file of randomly generated germline sequences.
#' @export
#'
#' @examples
#' data("ighv_hum_df")
#' data("ighd_hum_df")
#' data("ighj_hum_df")
#' germline_fasta = germline2fas(10)
germline2fas = function(family_num){
  indV = sample(x = 1:nrow(ighv_hum_df), family_num, replace = TRUE)
  indD = sample(x = 1:nrow(ighd_hum_df), family_num, replace = TRUE)
  indJ = sample(x = 1:nrow(ighj_hum_df), family_num, replace = TRUE)
  vseq = toupper(as.character(ighv_hum_df[[2]][indV]))
  dseq = toupper(as.character(ighd_hum_df[[2]][indD]))
  jseq = toupper(as.character(ighj_hum_df[[2]][indJ]))
  germline = paste(vseq, dseq, jseq, sep="")
  ids = c(1:family_num)
  fasta.df = data.frame(ids, germline)
  fasta = df2fas(fasta.df)
  return(fasta)
}

# germline_data = read.table('/data2/KMine/Simulation/germline.tsv', header=T, sep="\t")

#
# # Insertion
# ins = function(dna, loc){
#   left_dna = substring(dna, 1, loc)
#   len = nchar(dna)
#   right_dna = substring(dna, loc+1, len)
#   ran_ins_dna = sample(x=c("A","T","G","C"), 1, replace=TRUE, c(rep(1/4,4)))
#   ins_dna = paste(left_dna, ran_ins_dna, right_dna, sep = "")
#   return(ins_dna)
# }
#
# # Deletion
# del = function(dna, loc){
#   left_dna = substring(dna, 1, loc-1)
#   len = nchar(dna)
#   right_dna = substring(dna, loc+1, len)
#   del_dna = paste(left_dna, right_dna, sep = "")
#   return(del_dna)
# }
#
# # Substitution
# sub = function(dna, loc){
#   holding_char = substr(dna, loc, loc)
#   if(holding_char == "A")
#     substr(dna, loc, loc) = sample(x=c("T","G","C"), 1, replace=TRUE, c(15,70,15))
#   else if(holding_char=="T")
#     substr(dna, loc, loc) = sample(x=c("A","G","C"), 1, replace=TRUE, c(15,15,70))
#   else if(holding_char=="G")
#     substr(dna, loc, loc) = sample(x=c("A","T","C"), 1, replace=TRUE, c(70,15,15))
#   else if(holding_char=="C")
#     substr(dna, loc, loc) = sample(x=c("A","G","T"), 1, replace=TRUE, c(15,15,70))
#   return(dna)
# }
#
# # 每轮模拟过程
# each_round = function(input, mut = mut_ratio, not_mut = not_mut_ratio, mod = mode, start = cdr3_start, end = cdr3_end){
#   mut_n = length(input)
#   mut_dna = c()
#   for(kk in 1:mut_n){
#     cdr3_dna = input[kk]
#     act_n = sample(1:3, 1)
#     act_dna = rep(cdr3_dna, act_n)
#     for (ii in 1:act_n){ # 每轮中的每条序列
#       tmp_dna = act_dna[ii]
#       mut_loc = c()
#       mut_form = c()
#       for(cc in start:end){
#         CDR_mut = sample(x = c(0,1), 1, replace = TRUE, c(not_mut, mut))
#         if(CDR_mut == 1){
#           form = mod[sample(1:509, 1)]
#           mut_loc = c(mut_loc, cc)
#           mut_form = c(mut_form, form)
#         }
#       }
#       if(length(mut_loc) == 0)
#         mut_dna = c(mut_dna, tmp_dna)
#       else{
#         for(ll in 1:length(mut_loc)){
#           tmp.form = mut_form[ll]
#           tmp.loc = mut_loc[ll]
#           if(tmp.form == "ins"){
#             tmp_dna = ins(tmp_dna, tmp.loc)
#           }
#           else if(tmp.form == "del"){
#             tmp_dna = del(tmp_dna, tmp.loc)
#           }
#           else{
#             tmp_dna = sub(tmp_dna, tmp.loc)
#           }
#         }
#         mut_dna = c(mut_dna, tmp_dna)
#       }
#     }
#   }
#   return(mut_dna)
# }
#
#
# # 设置基本模拟参数
# load('/data2/KMine/Simulation/germline_pre.Rdata')
# n_true = 132
# data = germline_data[sample(1:nrow(germline_data), n_true, replace=F), ]
# mut_ratio = 0.01
# mode = c(rep("ins", 4), rep("del", 5), rep("sub", 500)) # 0.8:1:100
#
# filename = "~/changeo/132_0.01.fasta"
# annoname = "/data2/KMine/Simulation/132_0.001_igblast_db-pass_parse-select.tsv"
#
# final_dna = NULL
# for(mm in 1:nrow(data)){
#   # 第一轮
#   not_mut_ratio = 1 - mut_ratio
#   ori_dna = data[mm, "sequence"]
#   cdr3_start = data[mm, "v_sequence_end"] - 15
#   cdr3_end = data[mm, "j_sequence_start"] + 15
#   cdr3_len = cdr3_end - cdr3_start
#   act_n = 3
#   act_dna = rep(ori_dna, act_n)
#   mut_dna1 = c()
#   for (ii in 1:act_n){ # 每轮中的每条序列
#     tmp_dna = act_dna[ii]
#     mut_loc = c()
#     mut_form = c()
#     for(cc in cdr3_start:cdr3_end){
#       CDR_mut = sample(x = c(0,1), 1, replace = TRUE, c(not_mut_ratio, mut_ratio))
#       if(CDR_mut == 1){
#         form = mode[sample(1:509, 1)]
#         mut_loc = c(mut_loc, cc)
#         mut_form = c(mut_form, form)
#       }
#     }
#     if(length(mut_loc) == 0)
#       mut_dna1 = c(mut_dna1, tmp_dna)
#     else{
#       for(ll in 1:length(mut_loc)){
#         tmp.form = mut_form[ll]
#         tmp.loc = mut_loc[ll]
#         if(tmp.form == "ins"){
#           tmp_dna = ins(tmp_dna, tmp.loc)
#         }
#         else if(tmp.form == "del"){
#           tmp_dna = del(tmp_dna, tmp.loc)
#         }
#         else{
#           tmp_dna = sub(tmp_dna, tmp.loc)
#         }
#       }
#       mut_dna1 = c(mut_dna1, tmp_dna)
#     }
#   }
#
#   # 第二-七轮
#   mut_dna2 = each_round(mut_dna1)
#   mut_dna3 = each_round(mut_dna2)
#   mut_dna4 = each_round(mut_dna3)
#   mut_dna5 = each_round(mut_dna4)
#   mut_dna6 = each_round(mut_dna5)
#   mut_dna7 = each_round(mut_dna6)
#
#   clu = c(ori_dna, mut_dna1, mut_dna2, mut_dna3, mut_dna4, mut_dna5, mut_dna6, mut_dna7)
#   final_dna[[mm]] = clu
# }
#
# # 构建fasta并注释
# ids = c()
# seqs = c()
# for(i in 1:n_true){
#   tmp.seqs = final_dna[[i]]
#   seqs = c(seqs, tmp.seqs)
#   for(j in 1:length(tmp.seqs)){
#     tmp.ids = paste(i, '_', j, sep = '')
#     ids = c(ids, tmp.ids)
#   }
# }
# fasta.df = data.frame(ids,seqs)
# df2fa = function(data, file1){
#   ids = data$ids
#   seqs = data$seqs
#   res = paste0(">",ids,"\n",seqs)
#   write.table(res,
#               file = file1,
#               row.names = F,
#               col.names = F,
#               quote = F)
# }
# df2fa(fasta.df, filename)
# # fasta = dataframe2fas(fasta.df, file = filename)
#
# # 载入注释好的数据
# anno_data = read.table(annoname, header=T, sep="\t")
# anno_data = data.pre(anno_data)
# # anno_data = anno_data[,c("sequence_id","sequence","rev_comp","productive",
# #                           "v_call","d_call","j_call","sequence_alignment",
# #                           "germline_alignment","junction","junction_aa")]
# spe.ids = anno_data$sequence_id
# spe.label = unlist(lapply(spe.ids, function(x) unlist(strsplit(x, '_'))[1]))
# anno_data$label = spe.label
# save(anno_data, file = '/data2/KMine/Simulation/132_0.002_Hilary.Rdata')
#
# # 加入噪声序列
# # load("/data2/KMine/covid_sample_500k.Rdata")
# # noise.num = 20000
# # noise_loc = sample(1:nrow(sampledata), noise.num)
# # noi_data = sampledata[noise_loc, ]
# # noi_data = noi_data[,c("sequence_id","sequence","rev_comp","productive",
# #                        "v_call","d_call","j_call","sequence_alignment",
# #                        "germline_alignment","junction","junction_aa")]
# # noi_data$label = rep('Noise', noise.num)
# # save(noi_data, file = "/data2/KMine/Simulation/noise.Rdata")
#
# # 选择模拟簇的数量并合并噪声
# load("/data2/KMine/Simulation/noise.Rdata")
# load('/data2/KMine/Simulation/132_0.005.Rdata')
# # 0.002
# # anno_data = anno_data[-which(anno_data$label == '88'),]
# # anno_data = anno_data[-which(anno_data$label == '103'),]
#
# # 0.005
# # anno_data = anno_data[-which(anno_data$label == '18'),]
# # anno_data = anno_data[-which(anno_data$label == '50'),]
# # anno_data = anno_data[-which(anno_data$label == '121'),]
# # anno_data = anno_data[-which(anno_data$label == '40'),]
#
# # 0.01
# anno_data = anno_data[-which(anno_data$label == '75'),]
# manuname = "/data2/KMine/Simulation/50_0.01.Rdata"
# n_true = 50
# filt.label = names(table(anno_data$label))[intersect(which(as.numeric(table(anno_data$label)) < 80), which(as.numeric(table(anno_data$label)) > 30))]
# anno_data = anno_data[which(anno_data$label %in% filt.label), ]
# choose.label = sample(names(table(anno_data$label)), n_true, replace = F)
# choose_data = anno_data[which(anno_data$label %in% choose.label), ]
# manu_data = rbind(choose_data, noi_data)
# save(manu_data, file = manuname)
