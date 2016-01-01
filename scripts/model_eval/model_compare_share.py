
noparam_sections = ["subs_ml_model","subs_ml_model_options"]


def extract_model_free_param(modelfile) :
  """
  read modelfile and return param[param number] = (param_label,param_value)
  """

  print("warning -- temporary hack to skip cat_manager parameters -- check parameter sync")
  param_no=-1
  param_map={}

  for line in open(modelfile) :
    word = line.strip().split()

    # rougly determine if this is a parameter:
    lw = len(word)
    if lw < 3 : continue
    if lw > 5 : continue
    section=word[0]
    if section in noparam_sections : continue 

    train_bit=word[-2]
    if not ( train_bit == "0" or train_bit == "1" ) : continue
    is_train=bool(int(train_bit))

    label=word[1]
    param_val=word[-1]
    if is_train and label.find("edge_strength") == -1 :
      param_map[param_no] = (label,float(param_val))

    param_no += 1

  return param_map

