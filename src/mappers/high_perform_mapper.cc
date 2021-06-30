

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include "tradequeue_mapper.h"
#include "mappers/null_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

//TODO List: Things we can do.
//select initial processor in a smart way based on CELL -> PROCESSOR mapping (reduce number of instances) & Better NUMA/cache behaviour?
//Mapping tags? Inside Regent code need to use my_task:set_mapping_tag_id (and avoid default Mapper tags) but can let us do special things.
//Control Replication - will need to think about sharding factor. This is not yet done at all.


//What does the default mapper do. Write it here.
//If its an inner task map everything virtually if allowed. We will just let the default mapper deal with inner tasks.
//Otherwise, it looks to see if we've Cached the task result, and if so we reuse that mapping. We can do this always now instead of only for non-reduction tasks.
//Its safe entirely since the cache looks for a previous instance of the task on the same processor & task ID with the same task Hash.
//If that mapping doesn't exist or we can't reuse it for some reason, then we make instances for any regions which don't currently have space. The default Mapper
//always allocates memory for reductions, however we know we can just reuse them. This is something we need to be careful for - we only want to reuse the reduction
//instance if it "belongs" to the tasks' target processor.

////Mapper code here
class high_perform_mapper : public NullMapper
{
    public:
        high_perform_mapper(MapperRuntime *rt, Machine machine, Processor local,
                              const char *mapper_name);
        ~high_perform_mapper(void);
        virtual Mapper::MapperSyncModel get_mapper_sync_model(void) const;
        virtual void select_steal_targets(const MapperContext         ctx,
                                          const SelectStealingInput&  input,
                                                SelectStealingOutput& output);

    public:
//Task mapping calls
        virtual void select_tasks_to_map(const MapperContext          ctx,
                                         const SelectMappingInput&    input,
                                               SelectMappingOutput&   output);
        virtual void select_task_options(const MapperContext    ctx,
                                         const Task&            task,
                                               TaskOptions&     output);
        virtual void map_task(const MapperContext      ctx,
                              const Task&              task,
                              const MapTaskInput&      input,
                                    MapTaskOutput&     output);
        virtual void map_inline(const Mapping::MapperContext ctx,
                               const InlineMapping& inline_op,
                               const MapInlineInput& input, 
                               MapInlineOutput& output);
        virtual void slice_task(const Mapping::MapperContext ctx,
                                const Task& task,
                                const SliceTaskInput& input,
                                SliceTaskOutput& output);

    public:
        virtual void configure_context(const Mapping::MapperContext ctx,
                                       const Task&                  task, 
                                             ContextConfigOutput&   output);
        virtual void select_partition_projection(const MapperContext  ctx,
                                                 const Partition& partition,
                                                 const SelectPartitionProjectionInput& input,
                                                 SelectPartitionProjectionOutput& output);
        virtual void map_partition(const Mapping::MapperContext ctx,
                                   const Partition& partition,
                                   const MapPartitionInput& input,
                                         MapPartitionOutput& output);
        virtual void memoize_operation(const MapperContext  ctx,
                                       const Mappable&      mappable,
                                       const MemoizeInput&  input,
                                             MemoizeOutput& output);

        virtual void select_task_sources(const MapperContext        ctx,
                                         const Task&                task,
                                         const SelectTaskSrcInput&  input,
                                               SelectTaskSrcOutput& output);


    protected:
        virtual const char* get_mapper_name(void) const;
        virtual const char* create_mapper_name(void);
        std::map<std::tuple<Processor, IndexSpaceID, RegionTreeID>, std::vector<PhysicalInstance>> cell_reduction_instances;
        const char* mapper_name;
        virtual Processor select_initial_processor(const MapperContext ctx,
                                                   const Task&         task);

    //Only concerned with CPUs right now, other processor types can be added later.
    protected:
        //Who am I
        Processor local_proc;
        //System layout
        uint32_t total_nodes;
        AddressSpace node_id;
        std::vector<Processor> local_cpus;
        std::vector<Processor> remote_cpus;
        //For round-robining of tasks onto processor
        uint32_t next_local_cpu;
        Processor next_global_cpu;
        virtual Processor default_get_next_local_cpu();
        virtual Processor default_get_next_global_cpu();
        Machine::ProcessorQuery *global_cpu_query;

    protected:
        //Memory choice functionality
        virtual Memory get_best_memory(const Processor& target_proc);

    protected:
        //Maps from Unique cell to Processor.
        std::map<std::pair<IndexSpaceID, RegionTreeID>, Processor> cell_to_processor;

    protected:
        //Choose Physical Instances
        virtual PhysicalInstance choose_instance_default(const MapperContext ctx,
                                                         const RegionRequirement& req,
                                                         const Processor& target_proc);
        virtual void default_policy_select_constraints(MapperContext ctx,
                                                       LayoutConstraintSet &constraints,
                                                       Memory target_memory,
                                                       const RegionRequirement &req);
        virtual void select_soa_constraints(MapperContext ctx,
                                            LayoutConstraintSet &constraints,
                                            Memory target_memory,
                                            const RegionRequirement &req);
    protected:
         static inline bool physical_sort_func(
                         const std::pair<PhysicalInstance,unsigned> &left,
                         const std::pair<PhysicalInstance,unsigned> &right)
        { return (left.second < right.second); }
};

const char* high_perform_mapper::create_mapper_name(void)
{
    const char *result = (char*)malloc(128 * sizeof(char));
    snprintf(const_cast<char*>(result), 127, "High Performance Mapper for RegentParticleDSL");
    return result;
}


void high_perform_mapper::select_soa_constraints(MapperContext ctx,
                                                 LayoutConstraintSet &constraints,
                                                 Memory target_memory,
                                                 const RegionRequirement &req)
{
        // Our base default mapper will try to make instances of containing
        // all fields (in any order) laid out in SOA format to encourage
        // maximum re-use by any tasks which use subsets of the fields
//        constraints.add_constraint(SpecializedConstraint())
//          .add_constraint(MemoryConstraint(target_memory.kind()));
//
//        if (constraints.field_constraint.field_set.size() == 0)
//        {
//          // Normal instance creation
//          std::vector<FieldID> fields;
//          default_policy_select_constraint_fields(ctx, req, fields);
//          constraints.add_constraint(FieldConstraint(fields,false/*contiguous*/,
//                                                     false/*inorder*/));
//        }
//        if (constraints.ordering_constraint.ordering.size() == 0)
//        {
//          IndexSpace is = req.region.get_index_space();
//          Domain domain = runtime->get_index_space_domain(ctx, is);
//          int dim = domain.get_dim();
//          std::vector<DimensionKind> dimension_ordering(dim + 1);
//          for (int i = 0; i < dim; ++i)
//            dimension_ordering[i] =
//              static_cast<DimensionKind>(static_cast<int>(LEGION_DIM_X) + i);
//          dimension_ordering[dim] = LEGION_DIM_F;
//          constraints.add_constraint(OrderingConstraint(dimension_ordering,
//                                                        false/*contigous*/));
//        }
}

void high_perform_mapper::default_policy_select_constraints(MapperContext ctx,
                                                            LayoutConstraintSet &constraints,
                                                            Memory target_memory,
                                                            const RegionRequirement &req)
{
    //Use SoA layout for now as DefaultMapper creates
    select_soa_constraints(ctx, constraints, target_memory, req);
}

//TODO: Implement
high_perform_mapper::high_perform_mapper(MapperRuntime *rt, Machine machine, Processor local,
                             const char *mapper_name/*,
                             std::vector<Processor>* _procs_list,
                             std::vector<Memory>* _sysmems_list,
                             std::map<Memory, std::vector<Processor> >* _sysmem_local_procs,
                             std::map<Processor, Memory>* _proc_sysmems,
                             std::map<Processor, Memory>* _proc_regmems*/)
  : NullMapper(rt, machine),
    mapper_name(create_mapper_name()),
    local_proc(local),
    // procs_list(*_procs_list),
    // sysmems_list(*_sysmems_list),
    // sysmem_local_procs(*_sysmem_local_procs),
    // proc_sysmems(*_proc_sysmems)// ,
    // proc_regmems(*_proc_regmems)
    node_id(local.address_space()),
    next_local_cpu(0),
    global_cpu_query(NULL)
{
    //Setup processor map - taken from defaultmapper
    // Get all the processors and gpus on the local node
    Machine::ProcessorQuery all_procs(machine);
   for (Machine::ProcessorQuery::iterator it = all_procs.begin(); it != all_procs.end(); it++)
   {
        AddressSpace node = it->address_space();
        //Build local processor list first.
        if(node == node_id)
        {
            if(it->kind() == Processor::LOC_PROC)
            {
                local_cpus.push_back(*it);
            }
        }
        //Build global processor lists next.
        if(it->kind() == Processor::LOC_PROC)
        {
            if (node >= remote_cpus.size()){
                remote_cpus.resize(node+1, Processor::NO_PROC);
            }
            if(!remote_cpus[node].exists())
            {
                remote_cpus[node] = *it;
            }
        }
        
   }
   total_nodes = remote_cpus.size();
}

high_perform_mapper::~high_perform_mapper(void)
{
    free(const_cast<char*>(mapper_name));
}

const char* high_perform_mapper::get_mapper_name(void) const{
    return mapper_name;
}

//--------------------------------------------------------------------------
Mapper::MapperSyncModel high_perform_mapper::get_mapper_sync_model(void) const
    //--------------------------------------------------------------------------
{
  // Default mapper operates with the serialized re-entrant sync model
  return SERIALIZED_REENTRANT_MAPPER_MODEL;
}

//Select next CPUs, taken from Default Mapper
//--------------------------------------------------------------------------
Processor high_perform_mapper::default_get_next_local_cpu(void)
//--------------------------------------------------------------------------
{
    Processor result = local_cpus[next_local_cpu++];
    if (next_local_cpu == local_cpus.size())
        next_local_cpu = 0;
    return result;
}

//--------------------------------------------------------------------------
Processor high_perform_mapper::default_get_next_global_cpu(void)
//--------------------------------------------------------------------------
{
    if (total_nodes == 1)
        return default_get_next_local_cpu();
    if (!next_global_cpu.exists())
    {
        global_cpu_query = new Machine::ProcessorQuery(machine);
        global_cpu_query->only_kind(Processor::LOC_PROC);
        next_global_cpu = global_cpu_query->first();
    }
    Processor result = next_global_cpu;
    next_global_cpu = global_cpu_query->next(result);
    if (!next_global_cpu.exists())
    {
        delete global_cpu_query;
        global_cpu_query = NULL;
    }
  return result;
}

Processor high_perform_mapper::select_initial_processor(const MapperContext ctx,
                                                        const Task&         task)
{
    const char* task_name = task.get_task_name();
    if (strcmp(task_name, "pairwise_task") == 0 || strcmp(task_name, "self_task" ) == 0 || strcmp(task_name, "asym_pairwise_task") == 0)
    {
        //All these task launches should be index launched.
        assert(task.is_index_space);
        //Construct the unique identifier for the cell.
        std::pair<IndexSpaceID, RegionTreeID> uniq_id = std::pair<IndexSpaceID, RegionTreeID>{task.regions[0].region.get_index_space().get_id(),
                                                                                              task.regions[0].region.get_tree_id()};

        //Check if this has a processor chosen.
        std::map<std::pair<IndexSpaceID, RegionTreeID>, Processor>::const_iterator resulting_processor = cell_to_processor.find(uniq_id);
        if(resulting_processor != cell_to_processor.end())
        {
            return resulting_processor->second;
        }else
        {
            //Pick next processor globally
            Processor chosen_processor = default_get_next_global_cpu();
            cell_to_processor[uniq_id] = chosen_processor;
            return chosen_processor;
        }
    }

    //For all other task types just pick a processor
    if( task.is_index_space )
    {
        //For index space tasks, choose myself
        return local_proc;
    }

    //Otherwise, check depth a la DefaultMapper
    const int depth = task.get_depth();
    switch(depth)
    {
        case 0:
            //Top-level task. Keep this here
            return local_proc;
        
       case 1:
            //Distribute evenly around the machine
            return default_get_next_global_cpu();

        default:
            //Keep locally otherwise, assume at depth 1 these were distributed
            return default_get_next_local_cpu();

    }

    //We should never get here
    assert(false);
    return Processor::NO_PROC;
}

void high_perform_mapper::select_tasks_to_map(const MapperContext          ctx,
                                              const SelectMappingInput&    input,
                                                    SelectMappingOutput&   output)
{
    //Copy gather_perf mapper and just map everything!
    output.map_tasks.insert(input.ready_tasks.begin(),
                            input.ready_tasks.end());
}

void high_perform_mapper::configure_context(const Mapping::MapperContext ctx,
                                            const Task&                  task, 
                                                  ContextConfigOutput&   output)
{
    // defaults are fine?
    // Need to look what this actually does at some point
}

Memory high_perform_mapper::get_best_memory(const Processor& target_proc)
{
        Machine::MemoryQuery visible_memories(machine);
        visible_memories.has_affinity_to(target_proc);
        Memory best_memory = Memory::NO_MEMORY;
        unsigned best_bandwidth = 0;
        std::vector<Machine::ProcessorMemoryAffinity> affinity(1);
        for (Machine::MemoryQuery::iterator it = visible_memories.begin();
            it != visible_memories.end(); it++)
        {
            affinity.clear();
            machine.get_proc_mem_affinity(affinity, target_proc, *it,
                                          false /*not just local affinities*/);
            assert(affinity.size() == 1);
            if (!best_memory.exists() || (affinity[0].bandwidth > best_bandwidth)) {
                best_memory = *it;
                best_bandwidth = affinity[0].bandwidth;
            }
        }
        //TODO: Cache best memory per processor
        return best_memory;
}

PhysicalInstance high_perform_mapper::choose_instance_default(const MapperContext ctx,
                                                              const RegionRequirement& req,
                                                              const Processor& target_proc)
{
        Memory best_memory = get_best_memory(target_proc);
        LayoutConstraintSet constraints;
        FieldConstraint fc(req.privilege_fields, false /*!contiguous*/);
        constraints.add_constraint(fc);
        if (req.privilege == LEGION_REDUCE)
        {
            constraints.add_constraint(SpecializedConstraint(
                            LEGION_AFFINE_REDUCTION_SPECIALIZE, req.redop));
        }/*else{
            FieldConstraint fc(false );
            runtime->get_field_space_fields(ctx, req.region.get_field_space(),
                      fc.field_set);
            constraints.add_constraint(fc);
        }*/
                 
        std::vector<LogicalRegion> regions(1, req.region);
        Mapping::PhysicalInstance result;
        bool created;
        bool ok = runtime->find_or_create_physical_instance(ctx,
                                best_memory,
                                constraints,
                                regions,
                                result,
                                created);
        assert(ok);
        return result;
}

//Time to map task
void high_perform_mapper::map_task(const MapperContext      ctx,
                                   const Task&              task,
                                   const MapTaskInput&      input,
                                         MapTaskOutput&     output)
{
    //No postmapping of tasks right now
    output.postmap_task = false;

    //For now just run the task on the target processor
    output.target_procs.push_back(task.target_proc);

    const char* task_name = task.get_task_name();
    //Check the variant is ok (should be)
    std::vector<VariantID> valid_variants;
    runtime->find_valid_variants(ctx, task.task_id, valid_variants, task.target_proc.kind());
    assert(!valid_variants.empty());
    output.chosen_variant = valid_variants[0];

    //TODO We could cache this later.
    //If its a special type of task, do something special
    if( (strcmp(task_name, "asym_pairwise_task") == 0) || (strcmp(task_name, "self_task") == 0))
    {
        std::vector<std::set<FieldID> > missing_fields(task.regions.size());
        for(size_t i = 0; i < task.regions.size(); i++)
        {
            if(task.regions[i].privilege == LEGION_REDUCE)
            {
                //Check if we already made the PhysicalInstance for this
                //Construct a triplet to uniquely represent any reduction region
                std::tuple<Processor, IndexSpaceID, RegionTreeID> triple = std::tuple<Processor, IndexSpaceID, RegionTreeID>{task.target_proc,
                                                                            task.regions[i].region.get_index_space().get_id(),
                                                                            task.regions[i].region.get_tree_id()};
                //Create an iterator over the list of valid options
                std::map<std::tuple<Processor, IndexSpaceID, RegionTreeID>, std::vector<PhysicalInstance>>::const_iterator
                            find_reduc_instance = cell_reduction_instances.find(triple);
                //If we find a valid option
                if (find_reduc_instance != cell_reduction_instances.end())
                {
                    //Find it. This should always be of size 1 but we don't check
                    std::vector<PhysicalInstance> valid_instances;
                    for(std::vector<PhysicalInstance>::const_iterator it = find_reduc_instance->second.begin();
                            it != find_reduc_instance->second.end(); it++)
                    {
                        valid_instances.push_back(*it);
                    }

                    std::set<FieldID> valid_missing_fields;
                    runtime->filter_instances(ctx, task, i, output.chosen_variant,
                                              valid_instances, valid_missing_fields);
#ifndef NDEBUG
                    bool check =
#endif
                    runtime->acquire_and_filter_instances(ctx, valid_instances);
                    assert(check);

                    output.chosen_instances[i] = valid_instances;
                    missing_fields[i] = valid_missing_fields;
                    //If we got all the fields (which I think we should have) then we go to the next region
                    if (missing_fields[i].empty()){
                        continue;
                    }
                }else
                {

                    Mapping::PhysicalInstance inst;
                    RegionRequirement req = task.regions[i];
                    if (req.privilege_fields.size() != 0){
                        inst = choose_instance_default(ctx, req, task.target_proc);
                        runtime->set_garbage_collection_priority(ctx, inst, LEGION_GC_NEVER_PRIORITY);
                        output.chosen_instances[i].push_back(inst);
                        cell_reduction_instances[triple] = output.chosen_instances[i];
                    }
                    continue;
                }
                //We should never get here as we assume we find all the fields or we find none.
                assert(false);
            }else
            {
                Mapping::PhysicalInstance inst;
                RegionRequirement req = task.regions[i];

                if (req.privilege_fields.size() != 0){
                    inst = choose_instance_default(ctx, req, task.target_proc);
                    runtime->set_garbage_collection_priority(ctx, inst, -1);
                    output.chosen_instances[i].push_back(inst);
                }
            }
        }
    }else //For other types of task do the most naive option
    {
        for(size_t i = 0; i < task.regions.size(); i++)
        {
            Mapping::PhysicalInstance inst;
            RegionRequirement req = task.regions[i];
    
            if (req.privilege_fields.size() != 0){
                inst = choose_instance_default(ctx, req, task.target_proc);
                output.chosen_instances[i].push_back(inst);
            }
        }
    }
//    assert(false);
}

void high_perform_mapper::map_inline(const Mapping::MapperContext ctx,
                                     const InlineMapping& inline_op,
                                     const MapInlineInput& input, 
                                           MapInlineOutput& output)
{
    Mapping::PhysicalInstance inst;
    inst = choose_instance_default(ctx, inline_op.requirement, inline_op.parent_task->current_proc);
    output.chosen_instances.push_back(inst);
}


void high_perform_mapper::select_partition_projection(const Mapping::MapperContext  ctx,
                                                      const Partition& partition,
                                                      const SelectPartitionProjectionInput& input,
                                                      SelectPartitionProjectionOutput& output)
{
    // If we add repartitioning in we might want to redo this
    // If we have a complete partition then use it
    if (!input.open_complete_partitions.empty())
        output.chosen_partition = input.open_complete_partitions[0];
    else
        output.chosen_partition = LogicalPartition::NO_PART;
}


//Taken from Legions default_mapper right now
void high_perform_mapper::map_partition(const Mapping::MapperContext ctx,
                                        const Partition& partition,
                                        const MapPartitionInput& input,
                                              MapPartitionOutput& output)
{

    //   Mapping::PhysicalInstance inst;
//    inst = choose_instance(ctx, index_point, num_points,
//               partition.requirement);
      //Use the valid instances we were given if they still exist
      output.chosen_instances = input.valid_instances;
      if (!output.chosen_instances.empty())
        runtime->acquire_and_filter_instances(ctx,
                                          output.chosen_instances);
      // Now see if we have any fields which we still make space for
      std::vector<unsigned> to_erase;
      std::set<FieldID> missing_fields =
      partition.requirement.privilege_fields;
      for (std::vector<PhysicalInstance>::const_iterator it =
            output.chosen_instances.begin(); it !=
            output.chosen_instances.end(); it++)
      {
        if (it->get_location().kind() == Memory::GPU_FB_MEM) {
          // These instances are not supported yet (see Legion issue #516)
          to_erase.push_back(it - output.chosen_instances.begin());
        } else {
          it->remove_space_fields(missing_fields);
          if (missing_fields.empty())
            break;
        }
      }
      // Erase undesired instances
      for (std::vector<unsigned>::const_reverse_iterator it =
            to_erase.rbegin(); it != to_erase.rend(); it++)
      {
        output.chosen_instances.erase((*it) + output.chosen_instances.begin());
      }
      // If we've satisfied all our fields, then we are done
      if (missing_fields.empty())
        return;

      //Make an instance for our missing fields
       Memory best_memory = get_best_memory(partition.parent_task->current_proc);
      LayoutConstraintSet constraints;
      FieldConstraint fc(missing_fields, false /*!contiguous*/);
      constraints.add_constraint(fc);
      PhysicalInstance result;
      std::vector<LogicalRegion> regions(1, partition.requirement.region);
      bool created;
      bool ok = runtime->find_or_create_physical_instance(ctx,
                            best_memory,
                            constraints,
                            regions,
                            result,
                            created);
      assert(ok);
      output.chosen_instances.push_back(result);
}
//    // No constraints on mapping partitions
//      // Copy over all the valid instances, then try to do an acquire on them
//      // and see which instances are no longer valid
//      output.chosen_instances = input.valid_instances;
//      if (!output.chosen_instances.empty())
//        runtime->acquire_and_filter_instances(ctx,
//                                          output.chosen_instances);
//      // Now see if we have any fields which we still make space for
//      std::vector<unsigned> to_erase;
//      std::set<FieldID> missing_fields =
//        partition.requirement.privilege_fields;
//      for (std::vector<PhysicalInstance>::const_iterator it =
//            output.chosen_instances.begin(); it !=
//            output.chosen_instances.end(); it++)
//      {
//        if (it->get_location().kind() == Memory::GPU_FB_MEM) {
//          // These instances are not supported yet (see Legion issue #516)
//          to_erase.push_back(it - output.chosen_instances.begin());
//        } else {
//          it->remove_space_fields(missing_fields);
//          if (missing_fields.empty())
//            break;
//        }
//      }
//      // Erase undesired instances
//      for (std::vector<unsigned>::const_reverse_iterator it =
//            to_erase.rbegin(); it != to_erase.rend(); it++)
//        output.chosen_instances.erase((*it) + output.chosen_instances.begin());
//      // If we've satisfied all our fields, then we are done
//      if (missing_fields.empty())
//        return;
//      Machine::MemoryQuery visible_memories(machine);
//      Processor target_proc = partition.parent_task->current_proc;
//      visible_memories.has_affinity_to(target_proc);
//      Memory best_memory = Memory::NO_MEMORY;
//      unsigned best_bandwidth = 0;
//      std::vector<Machine::ProcessorMemoryAffinity> affinity(1);
//      for (Machine::MemoryQuery::iterator it = visible_memories.begin();
//          it != visible_memories.end(); it++)
//      {
//          affinity.clear();
//          machine.get_proc_mem_affinity(affinity, target_proc, *it,
//                                        false /*not just local affinities*/);
//          assert(affinity.size() == 1);
//          if (!best_memory.exists() || (affinity[0].bandwidth > best_bandwidth)) {
//              best_memory = *it;
//              best_bandwidth = affinity[0].bandwidth;
//          }
//      }
//            
//      //Not yet finished
//      Memory target_memory = best_memory;
//      //Select constraints
//      LayoutConstraintSet constraints;
//      default_policy_select_constraints(ctx, constraints, target_memory, req);
//      LayoutConstraintID our_layout_id =
//        runtime->register_layout(ctx, constraints);
//      //TODO: Could cache this again
//      LayoutConstraintSet creation_constraints =
//              runtime->find_layout_constraints(ctx, our_layout_id);
//      creation_constraints.add_constraint(FieldConstraint(missing_fields, false/*contig*/, false/*inorder*/));
//      output.chosen_instances.resize(output.chosen_instances.size()+1);
//
//      //We only have top level partitions for now I think
//      assert(!runtime->has_parent_logical_partition(ctx, output.chosen_instances.back()));
//      LogicalRegion target_region = partition.requirement;
//      bool tight_region_bounds = creation_constraints.specialized_constraint.is_exact()
//        || ((req.tag & DefaultMapper::EXACT_REGION) != 0);
//       std::vector<LogicalRegion> target_regions(1, target_region);
//       bool created = true;
//       if (force_new ||
//          ((req.privilege == LEGION_REDUCE) && (kind != COPY_MAPPING))) {
//        if (!runtime->create_physical_instance(ctx, target_memory,
//              creation_constraints, target_regions, output.chosen_instances.back(), true/*acquire*/,
//              0/*priority*/, tight_region_bounds, footprint))
//            {
//                assert(false); //Failed to allocate
//            }
//      } else {
//        if (!runtime->find_or_create_physical_instance(ctx,
//              target_memory, creation_constraints, target_regions, output.chosen_instances.back(), created,
//              true/*acquire*/, 0/*priority*/, tight_region_bounds, footprint))
//            {
//                assert(false); //Failed to allocate
//            }
//      }
//        if(created)
//        {
//            int priority = LEGION_GC_DEFAULT_PRIORITY;
//            runtime->set_garbage_collection_priority(ctx, output.chosen_instances.back(),priority);
//        }
//
//      assert(false);
//}

void high_perform_mapper::slice_task(const Mapping::MapperContext ctx,
          const Task& task, const SliceTaskInput& input,
          SliceTaskOutput& output)
{
    // even though we're going to map everything from the first processor,
    //  we need to slice apart all the points so that they can be mapped
    //  to different places
    for(Domain::DomainPointIterator dpi(input.domain); dpi; dpi.step())
      output.slices.push_back(TaskSlice(Domain(dpi.p, dpi.p),  default_get_next_global_cpu(),
                    false /*!recurse*/,
                    false /*!stealable*/));
}


void high_perform_mapper::select_task_options(const MapperContext    ctx,
                                              const Task&            task,
                                                    TaskOptions&     output)
{
    //Choose the initial processor according to the rules implemented in select_initial_processor 
    output.initial_proc = select_initial_processor(ctx, task); 

    //No stealing yet.
    output.stealable = false;
    //Map at target processor
    output.map_locally = false;
}


void high_perform_mapper::select_task_sources(const MapperContext        ctx,
                                              const Task&                task,
                                              const SelectTaskSrcInput&  input,
                                                   SelectTaskSrcOutput& output)
{
    //Copy the default mapper
      // For right now we'll rank instances by the bandwidth of the memory
      // they are in to the destination
      // TODO: consider layouts when ranking source  to help out the DMA system
      const PhysicalInstance target = input.target;
      const std::vector<PhysicalInstance> sources = input.source_instances;
      std::map<Memory,unsigned/*bandwidth*/> source_memories;
      Memory destination_memory = target.get_location();
      std::vector<MemoryMemoryAffinity> affinity(1);
      // fill in a vector of the sources with their bandwidths and sort them
      std::vector<std::pair<PhysicalInstance,
                          unsigned/*bandwidth*/> > band_ranking(sources.size());
      for (unsigned idx = 0; idx < sources.size(); idx++)
      {
        const PhysicalInstance &instance = sources[idx];
        Memory location = instance.get_location();
        std::map<Memory,unsigned>::const_iterator finder =
          source_memories.find(location);
        if (finder == source_memories.end())
        {
          affinity.clear();
          machine.get_mem_mem_affinity(affinity, location, destination_memory,
                       false /*not just local affinities*/);
          unsigned memory_bandwidth = 0;
          if (!affinity.empty()) {
            assert(affinity.size() == 1);
            memory_bandwidth = affinity[0].bandwidth;
          }
          source_memories[location] = memory_bandwidth;
          band_ranking[idx] =
            std::pair<PhysicalInstance,unsigned>(instance, memory_bandwidth);
        }
        else
          band_ranking[idx] =
            std::pair<PhysicalInstance,unsigned>(instance, finder->second);
      }
      // Sort them by bandwidth
      std::sort(band_ranking.begin(), band_ranking.end(), physical_sort_func);
      // Iterate from largest bandwidth to smallest
      for (std::vector<std::pair<PhysicalInstance,unsigned> >::
            const_reverse_iterator it = band_ranking.rbegin();
            it != band_ranking.rend(); it++)
        output.chosen_ranking.push_back(it->first);
}


void high_perform_mapper::memoize_operation(const Mapping::MapperContext ctx,
                                            const Mappable& mappable, 
                                            const MemoizeInput& input,
                                            MemoizeOutput& output)
{
    // memoize all the things
    output.memoize = true;
}



//Steal related functions.
void high_perform_mapper::select_steal_targets(const MapperContext         ctx,
                                               const SelectStealingInput&  input,
                                                SelectStealingOutput& output)
{
    //Not yet doing anything
}


static void create_mappers(Machine machine, Runtime *runtime, const std::set<Processor> &local_procs){
  for (std::set<Processor>::const_iterator it = local_procs.begin();
        it != local_procs.end(); it++)
  {
    high_perform_mapper* mapper = new high_perform_mapper(runtime->get_mapper_runtime(),
                                              machine, *it, "high_perform_mapper"/*,
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


//Next up TODO "select_tasks_to_map"
