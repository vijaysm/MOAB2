include ../make.common
INC_DIR = ../include
SRC_DIR = ../src
LIB_DIR = ../lib

VPATH = $(SRC_DIR):$(INC_DIR)

.phony: default shared

default: $(LIB_DIR)/libmhdf-pic.a
shared:  $(LIB_DIR)/libmhdf.so

.c.o: $<
	$(CC) $(PIC) $(DEBUG_FLAG) $(CLFAGS) -I$(INC_DIR) $(HDF5_INC) -o $@ -c $< 

$(LIB_DIR)/libmhdf-pic.a: $(OBJS)
	$(ARCHIVER) $@ $(OBJS)

$(LIB_DIR)/libmhdf.so: $(OBJS)
	$(SO_LINKER) $(OBJS) $(HDF5_LNK) -lm -o $@

clean:
	rm -f $(OBJS) $(LIB_DIR)/libmhdf-pic.a $(LIB_DIR)/libmhdf.so

