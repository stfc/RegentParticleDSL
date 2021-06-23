

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include "tradequeue_mapper.h"
#include "mappers/default_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

////Mapper code here
class tradequeueMapper : public DefaultMapper
{
    public:
    tradequeueMapper(MapperRuntime *rt, Machine machine, Processor local,
                const char *mapper_name);/*,
                std::vector<Processor>* procs_list,
                std::vector<Memory>* sysmems_list,
                std::map<Memory, std::vector<Processor> >* sysmem_local_procs,
                std::map<Processor, Memory>* proc_sysmems,
                std::map<Processor, Memory>* proc_regmems); TODO: Might want this stuff*/
    virtual Processor default_policy_select_initial_processor(
                                    MapperContext ctx, const Task &task);

    virtual void default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs);

      virtual void map_task(const MapperContext ctx,
                        const Task &task,
                        const MapTaskInput &input,
                        MapTaskOutput &output);
};

//TODO: Implement
//
tradequeueMapper::tradequeueMapper(MapperRuntime *rt, Machine machine, Processor local,
                             const char *mapper_name/*,
                             std::vector<Processor>* _procs_list,
                             std::vector<Memory>* _sysmems_list,
                             std::map<Memory, std::vector<Processor> >* _sysmem_local_procs,
                             std::map<Processor, Memory>* _proc_sysmems,
                             std::map<Processor, Memory>* _proc_regmems*/)
  : DefaultMapper(rt, machine, local, mapper_name) //,
    // procs_list(*_procs_list),
    // sysmems_list(*_sysmems_list),
    // sysmem_local_procs(*_sysmem_local_procs),
    // proc_sysmems(*_proc_sysmems)// ,
    // proc_regmems(*_proc_regmems)
{
}


//TODO: Implement.
Processor tradequeueMapper::default_policy_select_initial_processor( MapperContext ctx, const Task &task){
//    const char* task_name = task.get_task_name();
//    if (strcmp(task_name, "pairwise_task") == 0 || strcmp(task_name, "self_task") == 0){
//        for( unsigned idx=0; idx < task.regions.size(); idx++){
//            std::cout << task.regions[idx].region.get_index_space().get_id() << "\n";
//        }
  //      std::cout << "\n\n\n";
//    }
    return DefaultMapper::default_policy_select_initial_processor(ctx, task);
}

//TODO: Implement.
void tradequeueMapper::default_policy_select_target_processors(MapperContext ctx, const Task &task,std::vector<Processor> &target_procs){
    return DefaultMapper::default_policy_select_target_processors(ctx, task, target_procs);
}

//TODO: Implement.
void tradequeueMapper::map_task(const MapperContext ctx, const Task &task, const MapTaskInput &input, MapTaskOutput &output){

//    const char* task_name = task.get_task_name();
    //For pairwise or self tasks we map ourselves
//    if (strcmp(task_name, "pairwise_task") == 0 || strcmp(task_name, "self_task") == 0){
//        //Find the processor target kind and variant (as default Mapper)
//        Processor::Kind target_kind = task.target_proc.kind();
//        VariantInfo chosen = default_find_preferred_variant(task, ctx, true, true, target_kind);
//    }


    // For now for unknown/uncared task types we fall through to the default Mapper
    return DefaultMapper::map_task(ctx, task, input, output);
}






















static void create_mappers(Machine machine, Runtime *runtime, const std::set<Processor> &local_procs){
  for (std::set<Processor>::const_iterator it = local_procs.begin();
        it != local_procs.end(); it++)
  {
    tradequeueMapper* mapper = new tradequeueMapper(runtime->get_mapper_runtime(),
                                              machine, *it, "tradequeue_mapper"/*,
                                              procs_list,
                                              sysmems_list,
                                              sysmem_local_procs,
                                              proc_sysmems,
                                              proc_regmems*/);
    runtime->replace_default_mapper(mapper, *it);
  }
}






//TODO: Check this is right
void register_mappers()
{
  Runtime::add_registration_callback(create_mappers);
}
