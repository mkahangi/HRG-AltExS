CXX = g++
LD = g++
TARGET = EoSHRG

OBJDIR = ./build
SRCDIR = .
INCLUDEDIR = .

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
HEADERS = $(wildcard $(SRCDIR)/*.hpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
DEPFILES = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.d,$(SOURCES))

LIBS =
CXXFLAGS = -I$(INCLUDEDIR) -Wall -Wextra -g3 -O3 -std=gnu++20
LDFLAGS =  -g3
MAKEFLAGS += -j$(nproc)
RM = rm -f

.PHONY: all clean run redo build rerun format
.SECONDARY: $(DEPFILES)

all: $(TARGET)

$(TARGET): $(OBJECTS) | build
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJDIR)/%.d: $(SRCDIR)/%.cpp | build
	$(CXX) $(CXXFLAGS) -MM -MP -MT '$(OBJDIR)/$*.o' -MF $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/%.d | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p $(OBJDIR)

clean:
	$(RM) $(TARGET) $(OBJECTS) $(DEPFILES)

run:
	$(MAKE) build
	$(MAKE) $(TARGET)
	./$(TARGET)

redo:
	$(MAKE) clean
	$(MAKE) build
	$(MAKE)

rerun:
	$(MAKE) clean
	$(MAKE) build
	$(MAKE) run

format:
	clang-format -style llvm -i $(SOURCES) $(HEADERS)

# do not include dependency files for this recipes
ifneq ($(MAKECMDGOALS),clean, format)
-include $(DEPFILES)
endif
