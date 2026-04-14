
### Compile the IDL routine called 'create_goft_tables.pro' ==> then run the (similar) following command 
create_goft_tables, "fe_16_335.410"

### ======= If there are multiple close lines we can check it similar like the following ==========

## step 1: Run the following command
ch_synthetic, 93.92, 93.94, density=10^9, /goft, SNGL_ION='fe_18', output=ion

## Step 2: type the following in terminal 
help, ion

## we will see similar like this (which means there are two lines close to each other within the range of 93.92 and 93.94 Ang lines =======
# ** Structure <80109c08>, 17 tags, length=1952, data length=1912, refs=1:
#   LINES      STRUCT  -> <Anonymous> Array[2]
#   IONEQ_LOGT   FLOAT   Array[81]
#   IONEQ_NAME   STRING  'ch_adv_10-Jun-2025-16:56:49.ioneq'
#   IONEQ_REF    STRING  ''
#   WVL_LIMITS   FLOAT   Array[2]
#   MODEL_FILE   STRING  ' '
#   MODEL_NAME   STRING  'Constant density'
#   MODEL_NE    INT     -13824
#   MODEL_PE    FLOAT      0.00000
#   MODEL_TE    FLOAT      0.00000
#   WVL_UNITS    STRING  'Angstroms'
#   INT_UNITS    STRING  'erg cm+3 sr-1 s-1'
#   ADD_PROTONS   INT       1
#   DATE      STRING  'Tue Jun 10 16:57:12 2025'
#   VERSION     STRING  '11.0.2'
#   LOOKUP     INT       0
#   PHOTOEXCITATION INT       0

## Step 3: plot the multiple like the following in terminal 
plot, ion.ioneq_logt, ion.lines[0].goft , /ylog
oplot, ion.ioneq_logt, ion.lines[1].goft

## To see the wavelength info of the lines use the following command for the first line
print, ion.lines[0].wvl

### If there are multiple close lines, we can take a narrow wavelength limit to select a single line like below
create_goft_tables, "fe_18_93.932", wlim=0.0005
