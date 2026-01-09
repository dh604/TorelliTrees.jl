# Compute contributions with backups
function compute_contributions_backup(STs::Vector{StratumTree}, prefix::String; print_contributions::Bool=false)
  STB = ST_backup(STs, 0, prefix, print_contributions)
  _resume_compute_contributions_backup(STB) 
end

function load_and_resume_contributions_backup(filename::String)::Vector{StratumTree}
  STB = load_stratum_trees(filename)
  println("Backup file loaded successfully. Resuming computing contributions.")
  _resume_compute_contributions_backup(STB)
  return STB.STs
end

function _resume_compute_contributions_backup(STB::ST_backup)
  STs = STB.STs
  prefix = STB.prefix
  print_contributions = STB.print_contributions
  last_backup_time = now()
  n_trees = length(STs)
  for j in STB.last_index_with_cont_t+1:n_trees
    st = STs[j]
    cont = _compute_contribution(st)
    last_index_with_cont_t = j

    if print_contributions
      println("Contribution of tree number $j / , which is $st:")
      println(cont)
      println()
    end

    # create new backpup if necessary and delete old backup afterwards
    if now() - last_backup_time >= Day(1)
      
      _write_backup_and_delete_old(prefix, STs, print_contributions, last_index_with_cont_t)

      # reset timestamp for last backup.
      last_backup_time = now()
    end
  end

  # Also backup if everything is done.
  _write_backup_and_delete_old(prefix, STs, print_contributions, length(STs))
end

function _write_backup_and_delete_old(prefix::String, STs::Vector{StratumTree}, print_contributions::Bool, last_index_with_cont_t::Int64)
  backup_name_1 = prefix * "_A.jls"
  backup_name_2 = prefix * "_B.jls"
  sizediff = filesize(backup_name_2) - filesize(backup_name_1)
  backup_name = sizediff > 0 ? backup_name_1 : backup_name_2
  delete_name = sizediff > 0 ? backup_name_2 : backup_name_1

  STB = ST_backup(STs, last_index_with_cont_t, prefix, print_contributions)
  save_stratum_trees(backup_name, STB)
end