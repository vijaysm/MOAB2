include ../make.common
INC_DIR = ../include
SRC_DIR = ../src
LIB_DIR = ../lib

VPATH = $(SRC_DIR):$(INC_DIR)

.phony: default

default: $(LIB_DIR)/libmhdf.a

.c.o: $<
	$(CC) -DNDEBUG $(DEBUG_FLAG) $(CLFAGS) -I$(INC_DIR) $(HDF5_INC) -o $@ -c $< 

$(LIB_DIR)/libmhdf.a: $(OBJS)
	$(ARCHIVER) $@ $(OBJS)

clean:
	rm -f $(OBJS) $(LIB_DIR)/libmhdf.a
