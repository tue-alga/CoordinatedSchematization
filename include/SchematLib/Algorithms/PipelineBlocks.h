#ifndef SCHEMATLIB_ALGRITHMS_PIPELINEBLOCKS_H
#define SCHEMATLIB_ALGRITHMS_PIPELINEBLOCKS_H
#include <CppGeometryGUI/Utils/UIBuilder.h>

#include "SchematLib/MapSimplification/EdgeTrackingObserver.h"
#include "SchematLib/Models/BaseTypes.h"
#include "SchematLib/Models/EmbeddedGraph.h"

namespace SchematLib::Algorithms
{
    class PipelineBlock
    {
    public:
        virtual ~PipelineBlock() = default;
        // Observer to use for changes
        using SimpObs = SchematLib::MapSimplification::EdgeTrackingObserver;
        using Graph = Models::UndirectedEmbeddedGraph;
        virtual void apply(Graph& target, SimpObs& simplificationObserver, const Models::NT& bboxWidth, const Models::NT& bboxHeight) = 0;
        virtual void setupUi(CppGeometryGUI::UIBuilder& builder) = 0;
        virtual std::string type() const = 0;
        virtual std::string name() const { return m_name; };
        virtual void setName(const std::string& name)
        {
            m_name = name;
        };
        virtual std::string description() const = 0;
        virtual std::string label() const = 0;
    protected:
        std::string m_name;
    };

    template<typename PipelineBlockType>
    struct UiForPipelineBlock
    {
        void setupUi(CppGeometryGUI::UIBuilder& builder){}
    };

    // Block per pipeline element
}
#endif
